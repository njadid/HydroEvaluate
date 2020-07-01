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
from functools import partial
import itertools
import os

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
metric_path = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/metrics/Q'
metric_fn = 'Q_metrics.csv'

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
list2 = ['nRMSE', 'nMAE', 'timing', 'norm_bias', 'ppd', 'peak_qsim', 'qsim_vol', 'qobs_vol', 'pt_change_vol',
        'tim_peak']
header = ['year', 'sim_type', 'lid','agreementindex', 'bias', 'correlationcoefficient', 'covariance',
        'decomposed_mse', 'kge', 'log_p', 'lognashsutcliffe', 'mae', 'mse', 'nashsutcliffe', 'pbias',
        'rmse', 'rrmse', 'rsquared', 'rsr', 'volume_error'] + list2


def calc_metrics(args):
    lid,year = args

    metrics_result = []
    idx = np.where(lid_table['lid'] == lid)
    usgs_name = str(lid_table['USGS_id'][idx[0][0]])
    area = lid_table['Area'][idx[0][0]]

    # READ DATA TO A DICTIONARY 
    data = dict()
    # READ THE DATA FROM CSV FILES FOR EACH SIMULATION TYPE
    obs_fn = obs_path_format.format(lid=lid,year=year)
    try:
        data['obs'] = pd.read_csv(obs_fn, header=0, parse_dates=['dt'],
                           index_col=False)
        data['obs']['dt'] = pd.to_datetime(data['obs']['dt']).astype(int) / 10**9
    except:
        # print('err')
        return []
    avail_sim_types = []
    for sim_type in setup_list:
        sim_fn = sim_path_format.format(sim_type=sim_type,year=year,lid=lid)
        if not os.path.exists(sim_fn):
            continue
        avail_sim_types += [sim_type]
        data[sim_type] = pd.read_csv(sim_fn, header=0, parse_dates=['dt'],
                           index_col=False)
        data[sim_type]['dt'] = pd.to_datetime(data[sim_type]['dt']).astype(int) / 10**9
    # print(year,avail_sim_types)
    for sim_type in avail_sim_types:
        sub_data = pd.merge(data['obs'], data[sim_type], how='inner', on='dt', suffixes=('_obs', '_sim'))
        # print(sub_data)
        if len(sub_data) == 0:
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

        tim_peak = float(np.min(1.0 / 3600 * (t_sim_pk - t_obs_pk))*1e-9)  # in hours; +ve refers to delay
        RMSE = np.sqrt(np.mean((sub_data['Q_obs'] - sub_data['Q_sim']) ** 2))
        MAE = np.mean(np.abs(sub_data['Q_obs'] - sub_data['Q_sim']))
        nRMSE = RMSE / area  # Normalized by drainage area in square km
        nMAE = MAE / area  # Normalized by drainage area in square km

        metric_additional = [nRMSE, nMAE, timing, norm_bias, ppd, peak_qsim, qsim_vol, qobs_vol,
                                pt_change_vol, tim_peak]
        list_of_metrics = spotpy.objectivefunctions.calculate_all_functions(sub_data['Q_obs'],
                                                                            sub_data['Q_sim'])
        metrics_result += [[year, sim_type, lid] +[list_of_metrics[i][1] for i in range(len(list_of_metrics))] + metric_additional]
    return metrics_result


if __name__=='__main__':
    n_cpu = cpu_count()
    pool = ThreadPool(n_cpu)
    lid_list = lid_table['lid']
    permute_list = list(itertools.product(lid_list,year_list))
    arg_pool = [list(x) for x in permute_list]
    t_0 = datetime.now()
    results = pool.map(calc_metrics, arg_pool)
    t_1 = datetime.now()
    total_time = (t_1-t_0).total_seconds()/60
    print('Total Elapsed Time = ' + "%.2f" %total_time  + ' minutes.')
    results_all = [sublist_L2 for sublist_L1 in results for sublist_L2 in sublist_L1]
    metrics_df = pd.DataFrame(results_all, columns=header)
    metrics_df.dropna(inplace=True, axis=1)
    metrics_df = metrics_df.replace([np.inf, -np.inf], np.nan)
    metrics_df.dropna(inplace=True, axis=0)
    csv_fn = os.path.join(metric_path,metric_fn)
    if not os.path.exists(metric_path):
        os.makedirs(metric_path)
    metrics_df.to_csv(csv_fn, float_format='%.3f',index=False)
    print('Finished metric calculations!')