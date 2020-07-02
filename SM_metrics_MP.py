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
def xcorr(sim_new, obs, norm=True, maxlags=10):
    c = np.correlate(sim_new, obs, mode=2)
    Nx = np.array(sim_new).shape[0]
    if norm: c /= np.sqrt(np.dot(sim_new, sim_new) * np.dot(obs, obs))
    if maxlags is None: maxlags = Nx - 1
    lags = np.arange(-maxlags, maxlags + 1)
    c = c[Nx - 1 - maxlags:Nx + maxlags]
    return c, lags


# CFS TO CMS CONVERSION FACTOR
cfs2cms = 0.0283168

# COLUMN DELIMITER FOR USGS_GAGES
meas_delimiter = ','
dt_col = 'dt'
# DATA LOCATIONS FOR RESULTS (res_data_dir) AND GAGE OBSERVATION DATA (obs_data_dir)
res_data_dir = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/SM'
obs_data_dir = '/Users/njadidoleslam/data/sm_gages'

metric_path = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/metrics/SM'
metric_fn = 'SM_metrics.csv'

# READ THE LIST OF LINK_ID'S TO BE ANALYZED
lid_table = pd.DataFrame(pd.read_csv('sm_gage_list.csv'))

# lid = 74725

#List of simulation types 
setup_list = ['_10009_dist_p5', '_10009_dist_p50', '_10009_dist_p95',  '_10008_dist_p5', '_10008_dist_p50', '_10008_dist_p95', '_10006_dist_p5', '_10006_dist_p50', '_10006_dist_p95', '_10008_dist_matched_p5']


year_list = [2016, 2017, 2018, 2019]
depths = ['sm_5','sm_10','sm_20','sm_50']
sim_path_format = '/Users/njadidoleslam/hlm_dev/asynch_richards/results/SM/{sim_type}/{year}/{lid}.csv'
obs_path_format = '/Users/njadidoleslam/data/sm_gages/{year}/{lid}.csv'
# DEFINE A LIST OF METRICS
list2 = ['norm_bias', 'ref_std', 'sim_std']
header = ['year', 'sim_type', 'lid','agreementindex', 'bias', 'correlationcoefficient', 'covariance',
        'decomposed_mse', 'kge', 'log_p', 'lognashsutcliffe', 'mae', 'mse', 'nashsutcliffe', 'pbias',
        'rmse', 'rrmse', 'rsquared', 'rsr', 'volume_error'] + list2


def calc_metrics(args):
    lid,year = args
    metrics_result = []
    # READ DATA TO A DICTIONARY 
    data = dict()
    # READ THE DATA FROM CSV FILES FOR EACH SIMULATION TYPE
    obs_fn = obs_path_format.format(lid=lid,year=year)
    try:
        data['obs'] = pd.read_csv(obs_fn, header=0, parse_dates=[dt_col],
                           index_col=False)
        # data['obs'][dt_col] = pd.to_datetime(data['obs'][dt_col]).astype(int) / 10**9
    except:
        # print('err')
        return []
    avail_sim_types = []
    for i,sim_type in enumerate(setup_list):
        sim_fn = sim_path_format.format(sim_type=sim_type,year=year,lid=lid)
        if not os.path.exists(sim_fn):
            continue
        avail_sim_types += [sim_type]
        data[sim_type] = pd.read_csv(sim_fn, header=0, parse_dates=[dt_col],
                           index_col=False)
    for sim_type in avail_sim_types:

        sub_data = pd.merge(data['obs'], data[sim_type], how='inner', on=dt_col, suffixes=('_obs', '_sim'))
        sub_data.dropna(inplace=True)

        for depth in depths:

            sub_data.dropna(inplace=True)
            obs_name = depth + '_obs'
            sim_name = depth + '_sim'
            if sub_data is None or len(sub_data) == 0:
                continue

            ref_mean = np.mean(sub_data[obs_name])
            ref_std = np.std(sub_data[obs_name])
            sim_mean = np.mean(sub_data[sim_name])
            sim_std = np.std(sub_data[sim_name])
            tsim_unix = sub_data[dt_col]
            obs_vol = np.trapz(sub_data[obs_name], sub_data[dt_col])
            sim_vol = np.trapz(sub_data[sim_name], sub_data[dt_col])

            pt_change_vol = (sim_vol - obs_vol) * 100 / obs_vol

            norm_bias = (sim_mean - ref_mean) / ref_std
            peak_ref = np.max(sub_data[obs_name])
            t_obs_pk = tsim_unix[sub_data[obs_name] == peak_ref]
            peak_idx_sim = np.argmin(tsim_unix - t_obs_pk)
            peak_sim = np.max(sub_data.loc[peak_idx_sim - 5 * 24:peak_idx_sim + 5 * 24][sim_name])
            ppd = (peak_sim - peak_ref) / peak_ref
            t_sim_pk = tsim_unix[sub_data[sim_name] == peak_sim]
            tim_peak = float(np.min(1.0 / 3600 * (t_sim_pk - t_obs_pk))*1e-9)  # in hours; +ve refers to delay
            metric_additional = [norm_bias, ref_std, sim_std]
            list_of_metrics = spotpy.objectivefunctions.calculate_all_functions(sub_data[obs_name],
                                                                                sub_data[sim_name])
            metrics_result += [[year, sim_type+'_'+depth, lid] +[list_of_metrics[i][1] for i in range(len(list_of_metrics))] + metric_additional]
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
    total_time = (t_1-t_0)/60
    print('Total Elapsed Time = ' + "%.2f" % total_time.total_seconds() + ' seconds.')
    results_all = [sublist_L2 for sublist_L1 in results for sublist_L2 in sublist_L1]
    metrics_df = pd.DataFrame(results_all, columns=header)
    metrics_df.dropna(inplace=True, axis=1)
    metrics_df = metrics_df.replace([np.inf, -np.inf], np.nan)
    metrics_df.dropna(inplace=True, axis=0)
    csv_fn = os.path.join(metric_path,metric_fn)
    if not os.path.exists(metric_path):
        os.makedirs(metric_path)
    metrics_df.to_csv(csv_fn, float_format='%.3f',index=False)
    print('Finished SM metrics calculations!')
