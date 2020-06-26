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


# A FUNCTION THAT ACCOUNTS FOR THE DAYLIGHT SAVING TIMES IN THE USGS STREAMFLOW TIMESTAMPS
def cst_cdt_con(tz_cd):
    tz_off = []
    if tz_cd == "CDT":
        tz_off = '-5'
        tz_off = 5 * 3600
    elif tz_cd == "CST":
        # tz_off = '-6'
        tz_off = 6 * 3600
    return tz_off


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
meas_delimiter = '\t'

# DATA LOCATIONS FOR RESULTS (res_data_dir) AND USGS GAGE DATA (usgs_data_dir)
res_data_dir = 'E:/APPFiles/PostGreSql/EnterpriseDB-ApacheHTTPD/apache/www/smap/data/results'
usgs_data_dir = 'E:/APPFiles/PostGreSql/EnterpriseDB-ApacheHTTPD/apache/www/smap/data/usgs_gages_long'

# READ THE LIST OF LINK_ID'S TO BE ANALYZED
lid_table = pd.DataFrame(pd.read_csv('lid_usgs.csv'))

# lid = 74725

# DEFINE A LIST OF METRICS
list1 = ['lid', 'USGS_id', 'USGS_name', 'Area', 'Slope', 't_trav', 'lat', 'lon', 'ifis_id', 'num_up_gag']
list2 = ['nRMSE', 'nMAE', 'timing', 'norm_bias', 'ppd', 'peak_qsim', 'qsim_vol', 'qobs_vol', 'pt_change_vol',
         'tim_peak']
header = ['year', 'sim_type', 'agreementindex', 'bias', 'correlationcoefficient', 'covariance',
          'decomposed_mse', 'kge', 'log_p', 'lognashsutcliffe', 'mae', 'mse', 'nashsutcliffe', 'pbias',
          'rmse', 'rrmse', 'rsquared', 'rsr', 'volume_error'] + list1 + list2


def calc_metrics(lid):
    metrics_result = []
    idx = np.where(lid_table['lid'] == lid)
    usgs_name = str(lid_table['USGS_id'][idx[0][0]])
    lat = lid_table['Y'][idx[0][0]]
    lon = lid_table['X'][idx[0][0]]
    area = lid_table['Area'][idx[0][0]]

    # READ THE DATA FROM CSV FILES FOR EACH SIMULATION TYPE
    df_assim = pd.read_csv(path.join(res_data_dir, 'csv_assim/' + str(lid) + '.csv'), header=0, parse_dates=['dt'],
                           index_col=False)
    df_assim_norm = pd.read_csv(path.join(res_data_dir, 'csv_assim_norm/' + str(lid) + '.csv'), header=0,
                                parse_dates=['dt'], index_col=False)
    df_assim_norm_v2 = pd.read_csv(path.join(res_data_dir, 'csv_smap_l3_v2/' + str(lid) + '.csv'), header=0,
                                   parse_dates=['dt'], index_col=False)
    df_ol = pd.read_csv(path.join(res_data_dir, 'csv_ol/' + str(lid) + '.csv'), header=0, parse_dates=['dt'],
                        index_col=False)
    try:
        df_ol1 = pd.read_csv(path.join(res_data_dir, 'out_ol_new/lid_O_' + str(lid) + '.csv'), header=0,
                             parse_dates=['dt'], index_col=False)
        df_ol1 = df_ol1[['Q', 'dt']]
        df_ol = df_ol.append(df_ol1, sort=True, ignore_index=True)


    except:
        print('There is no data for this station for the simulation from 2010-2014')

    df_assim_smos = pd.read_csv(path.join(res_data_dir, 'csv_smos/' + str(lid) + '.csv'), header=0, parse_dates=['dt'])
    df_meas = pd.read_csv(path.join(usgs_data_dir, 'fifteen_minute', 'usgs_0' + usgs_name + '.txt'), usecols=[2, 3, 4],
                          names=['dt', 'tz', 'Q'], header=0, parse_dates=['dt'], delimiter=meas_delimiter)

    # CONVERT DATE TIME TO ONE TYPE OF DATETIME IN PYTHON
    df_assim['dt'] = [calendar.timegm(df_assim['dt'][i].timetuple()) for i in range(len(df_assim['dt']))]
    df_assim_norm['dt'] = [calendar.timegm(df_assim_norm['dt'][i].timetuple()) for i in range(len(df_assim_norm['dt']))]
    df_assim_norm_v2['dt'] = [calendar.timegm(df_assim_norm_v2['dt'][i].timetuple()) for i in
                              range(len(df_assim_norm_v2['dt']))]
    df_assim_smos['dt'] = [calendar.timegm(df_assim_smos['dt'][i].timetuple()) for i in range(len(df_assim_smos['dt']))]

    df_ol['dt'] = [calendar.timegm(df_ol['dt'][i].timetuple()) for i in range(len(df_ol['dt']))]
    df_ol.sort_values(by=['dt'], inplace=True)
    df_meas['dt'] = [calendar.timegm(df_meas['dt'][i].timetuple()) + cst_cdt_con(df_meas['tz'][i]) for i in
                     range(len(df_meas['dt']))]

    df_meas['Q'] = pd.to_numeric(df_meas['Q'], errors='coerce')
    df_meas.dropna(0, inplace=True)
    df_meas['Q'] *= cfs2cms

    data = dict()
    df_names = ['Q', 'Q_assim', 'Q_ol', 'Q_assim_norm', 'Q_assim_norm_v2', 'df_assim_smos']
    data['Q'] = df_meas
    data['Q_assim'] = df_assim
    data['Q_ol'] = df_ol
    data['Q_assim_norm'] = df_assim_norm
    data['Q_assim_norm_v2'] = df_assim_norm_v2
    data['Q_assim_smos'] = df_assim_smos

    for simulation in ['Q_ol', 'Q_assim', 'Q_assim_norm', 'Q_assim_norm_v2', 'Q_assim_smos']:
        result = pd.merge(data['Q'], data[simulation], how='inner', on='dt', suffixes=('', '_sim'))
        result['dt_datetime'] = [datetime.utcfromtimestamp(result['dt'][i]) for i in range(len(result['dt']))]
        result.columns = [u'dt', u'tz', u'Q_obs', u'Q_sim', 'dt_datetime']
        if not len(result) == 0:
            for year in [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]:
                subset_year = result[result['dt_datetime'].dt.year == year]
                if not len(subset_year) == 0:
                    qref_mean = np.mean(subset_year['Q_obs'])
                    qref_std = np.std(subset_year['Q_obs'])
                    qsim_mean = np.mean(subset_year['Q_sim'])
                    qsim_std = np.std(subset_year['Q_sim'])
                    tsim_unix = subset_year['dt']
                    qobs_vol = np.trapz(subset_year['Q_obs'], subset_year['dt'])
                    qsim_vol = np.trapz(subset_year['Q_sim'], subset_year['dt'])

                    pt_change_vol = (qsim_vol - qobs_vol) * 100 / qobs_vol

                    norm_bias = (qsim_mean - qref_mean) / qref_std
                    cor, lag = xcorr(subset_year['Q_sim'], subset_year['Q_obs'], maxlags=384)
                    cor_max = np.max(cor)
                    timing = round(lag[np.where(cor == cor_max)[0][0]] / 4)

                    peak_qref = np.max(subset_year['Q_obs'])
                    t_obs_pk = tsim_unix[subset_year['Q_obs'] == peak_qref]
                    peak_idx_sim = np.argmin(tsim_unix - t_obs_pk)
                    peak_qsim = np.max(subset_year.loc[peak_idx_sim - 5 * 24:peak_idx_sim + 5 * 24]['Q_sim'])
                    ppd = (peak_qsim - peak_qref) / peak_qref

                    t_sim_pk = tsim_unix[subset_year['Q_sim'] == peak_qsim]

                    tim_peak = np.min(1.0 / 3600 * (t_sim_pk - t_obs_pk))  # in hours; +ve refers to delay
                    RMSE = np.sqrt(np.mean((subset_year['Q_obs'] - subset_year['Q_sim']) ** 2))
                    MAE = np.mean(np.abs(subset_year['Q_obs'] - subset_year['Q_sim']))
                    nRMSE = RMSE / area  # Normalized by drainage area in square km
                    nMAE = MAE / area  # Normalized by drainage area in square km

                    metric_additional = [nRMSE, nMAE, timing, norm_bias, ppd, peak_qsim, qsim_vol, qobs_vol,
                                         pt_change_vol, tim_peak]
                    list_of_metrics = spotpy.objectivefunctions.calculate_all_functions(subset_year['Q_obs'],
                                                                                        subset_year['Q_sim'])
                    metrics_result += [[year, simulation] +
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
