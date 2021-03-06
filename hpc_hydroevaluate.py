# HYDROGRAPH EVALUATION AND COMPARISON WITH METRICS FOR SOIL MOISTURE ASSIMILATED RESULTS FROM HLM MODEL.
# THE SIMULATIONS ARE EVALUATED SEPARATELY FOR EACH YEAR, STATION AND SIMULATION TYPE

# IMPORT NECESSARY LIBRARIES todo: CAN BE LESS THAN THIS AND IT CAN BE OPTIMIZED
from os import path
import pandas as pd
import numpy as np
import calendar
from datetime import datetime
import spotpy
import json
# Import Multithreading Libraries
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count



# A SEPARATE FUNCTION FOR FINDING CROSS-CORRELATION AND LAG TIME
def xcorr(qsim_new, qobs, norm=True, maxlags=384):
    c = np.correlate(qsim_new, qobs, mode=2)
    Nx = np.array(qsim_new).shape[0]
    if norm: c /= np.sqrt(np.dot(qsim_new, qsim_new) * np.dot(qobs, qobs))
    if maxlags is None: maxlags = Nx - 1
    lags = np.arange(-maxlags, maxlags + 1)
    c = c[Nx - 1 - maxlags:Nx + maxlags]
    return c, lags


# CFS TO CMS CONVERSION FACTOR
cfs2cms = 0.0283168

# COLUMN DELIMITER FOR USGS_GAGES
meas_delimiter = ','

# DATA LOCATIONS FOR RESULTS (res_data_dir) AND USGS GAGE DATA (usgs_data_dir)
res_data_dir = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/Q'
usgs_data_dir = '/Users/njadidoleslam/data/usgs_obs'

# READ THE LIST OF LINK_ID'S TO BE ANALYZED
lid_table = pd.DataFrame(pd.read_csv('lid_usgs.csv'))

# lid = 74725

#List of simulation types 
setup_list = ['_254_mrms','_10009_dist_p5', '_10009_dist_p50', '_10009_dist_p95', 
            '_10008_dist_p5', '_10008_dist_p50', '_10008_dist_p95',
            '_10006_dist_p5', '_10006_dist_p50', '_10006_dist_p95', '_10008_dist_matched_p5']

#year_list = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]
year_list = [2016, 2017, 2018, 2019]
sim_path_format = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/Q/{sim_type}/{year}/{lid}.csv'
obs_path_format = '/Users/njadidoleslam/data/usgs_obs/{year}/{lid}.csv'
# DEFINE A LIST OF METRICS
list1 = ['lid', 'USGS_id', 'USGS_name', 'Area', 'Slope', 't_trav', 'lat', 'lon', 'ifis_id', 'num_up_gag']
list2 = ['nRMSE', 'nMAE', 'timing', 'norm_bias', 'ppd', 'peak_qsim', 'qsim_vol', 'qobs_vol', 'pt_change_vol',
         'tim_peak']
header = ['year', 'sim_type', 'agreementindex', 'bias', 'correlationcoefficient', 'covariance',
          'decomposed_mse', 'kge', 'log_p', 'lognashsutcliffe', 'mae', 'mse', 'nashsutcliffe', 'pbias',
          'rmse', 'rrmse', 'rsquared', 'rsr', 'volume_error'] + list1 + list2


def calc_metrics(lid,year):
    metrics_result = []
    idx = np.where(lid_table['lid'] == lid)
    usgs_name = str(lid_table['USGS_id'][idx[0][0]])
    area = lid_table['Area'][idx[0][0]]

    # READ DATA TO A DICTIONARY 
    data = dict()
    # READ THE DATA FROM CSV FILES FOR EACH SIMULATION TYPE
    obs_fn = obs_path_format.format(lid=lid,year=year)
    data['Q'] = pd.read_csv(obs_fn, header=0, parse_dates=['dt'],
                           index_col=False)
    for sim_type in setup_list:
        sim_fn = sim_path_format.format(sim_type=sim_type,year=year,lid=lid)
        data[sim_type] = pd.read_csv(sim_fn, header=0, parse_dates=['dt'],
                           index_col=False)

    for sim_type in setup_list:
        sub_data = pd.merge(data['Q'], data[sim_type], how='inner', on='dt', suffixes=('', '_sim'))
        sub_data['dt_datetime'] = [datetime.utcfromtimestamp(result['dt'][i]) for i in range(len(result['dt']))]
        sub_data.columns = [u'dt', u'Q_obs', u'Q_sim', 'dt_datetime']
        if not len(sub_data) == 0:
            continue
        qref_mean = np.mean(sub_data['Q_obs'])
        qref_std = np.std(sub_data['Q_obs'])
        qsim_mean = np.mean(sub_data['Q_sim'])
        qsim_std = np.std(sub_data['Q_sim'])
        tsim_unix = sub_data['dt']
        qobs_vol = np.trapz(sub_data['Q_obs'], sub_data['dt'])
        qsim_vol = np.trapz(sub_data['Q_sim'], sub_data['dt'])

        pt_change_vol = (qsim_vol - qobs_vol) * 100 / qobs_vol

        norm_bias = (qsim_mean - qref_mean) / qref_std
        cor, lag = xcorr(sub_data['Q_sim'], sub_data['Q_obs'], maxlags=384)
        cor_max = np.max(cor)
        timing = round(lag[np.where(cor == cor_max)[0][0]] / 4)

        peak_qref = np.max(sub_data['Q_obs'])
        t_obs_pk = tsim_unix[sub_data['Q_obs'] == peak_qref]
        peak_idx_sim = np.argmin(tsim_unix - t_obs_pk)
        peak_qsim = np.max(sub_data.loc[peak_idx_sim - 5 * 24:peak_idx_sim + 5 * 24]['Q_sim'])
        ppd = (peak_qsim - peak_qref) / peak_qref

        t_sim_pk = tsim_unix[sub_data['Q_sim'] == peak_qsim]

        tim_peak = np.min(1.0 / 3600 * (t_sim_pk - t_obs_pk))  # in hours; +ve refers to delay
        RMSE = np.sqrt(np.mean((sub_data['Q_obs'] - sub_data['Q_sim']) ** 2))
        MAE = np.mean(np.abs(sub_data['Q_obs'] - sub_data['Q_sim']))
        nRMSE = RMSE / area  # Normalized by drainage area in square km
        nMAE = MAE / area  # Normalized by drainage area in square km

        metric_additional = [nRMSE, nMAE, timing, norm_bias, ppd, peak_qsim, qsim_vol, qobs_vol,
                                pt_change_vol, tim_peak]
        list_of_metrics = spotpy.objectivefunctions.calculate_all_functions(sub_data['Q_obs'],
                                                                            sub_data['Q_sim'])
        metrics_result += [[year, sim_type] +
                            [list_of_metrics[i][1] for i in range(len(list_of_metrics))] +
                            [lid_table[x][idx[0][0]] for x in list1] + metric_additional]
    return metrics_result


# LOOP OVER THE GAGES INSIDE THE GAGE LIST: lid_table
# lid: LINK_ID OF THE MODEL
lid_list = lid_table['lid']
t_0 = datetime.now()
# Multithread computation starts
n_cpu = cpu_count()
pool = ThreadPool(n_cpu - 1)
results = pool.map(calc_metrics, lid_list)
# Multithread computation ends
t_1 = datetime.now()
total_time = (t_1-t_0)/60
print('Total Elapsed Time = ' + "%.2f" % total_time + ' Minutes.')
results_all = [sublist_L2 for sublist_L1 in results for sublist_L2 in sublist_L1]
df_final = pd.DataFrame(results_all, columns=header)
df_final1 = df_final
df_final1.dropna(inplace=True, axis=1)
df_final1 = df_final.replace([np.inf, -np.inf], np.nan)
df_final1.dropna(inplace=True, axis=0)
df_final1.to_csv('df_final_new.csv')

# metrics_dict = df_final1.to_dict()
# js = (json.dumps(metrics_dict).replace(u'<', u'\\u003c')
#       .replace(u'>', u'\\u003e')
#       .replace(u'&', u'\\u0026')
#       .replace(u"'", u'\\u0027'))
#
# fp = open('multithread_metrics.json', 'w')
# # write to json file
# fp.write(js)
# fp.close()
