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
import os
# Import Multithreading Libraries
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count

# CFS TO CMS CONVERSION FACTOR
cfs2cms = 0.0283168


def cst_cdt_con(tz_cd):
    tz_off = []
    if tz_cd == "CDT":
        tz_off = '-5'
        tz_off = 5 * 3600
    elif tz_cd == "CST":
        # tz_off = '-6'
        tz_off = 6 * 3600
    return tz_off


# COLUMN DELIMITER FOR USGS_GAGES
meas_delimiter = ','

# DATA LOCATIONS FOR RESULTS (res_data_dir) AND USGS GAGE DATA (usgs_data_dir)
res_data_dir = 'E:/APPFiles/PostGreSql/EnterpriseDB-ApacheHTTPD/apache/www/smap-new/data/results'
usgs_data_dir = 'E:/APPFiles/PostGreSql/EnterpriseDB-ApacheHTTPD/apache/www/smap/data/usgs_gages_long'
out_folder = 'E:\\platformdata'
# READ THE LIST OF LINK_ID'S TO BE ANALYZED
lid_table = pd.DataFrame(pd.read_csv('lid_usgs.csv'))


def calc_metrics(lid):
    metrics_result = []
    idx = np.where(lid_table['lid'] == lid)
    usgs_name = '0' + str(lid_table['USGS_id'][idx[0][0]])
    lat = lid_table['Y'][idx[0][0]]
    lon = lid_table['X'][idx[0][0]]
    area = lid_table['Area'][idx[0][0]]
    df_assim_norm_v2 = pd.read_csv(path.join(res_data_dir, 'csv_smap_l3_v2/' + str(lid) + '.csv'), header=0,
                                   parse_dates=['dt'], index_col=False)
    # Openloop data
    df_ol = pd.read_csv(path.join(res_data_dir, 'csv_ol/' + str(lid) + '.csv'), header=0, parse_dates=['dt'],
                        index_col=False)
    # USGS observations
    df_meas = pd.read_csv(path.join(usgs_data_dir, 'fifteen_minute1', 'usgs_' + usgs_name + '.txt'),
                          names=['dt', 'tz', 'Q'], header=0, parse_dates=['dt'], delimiter=meas_delimiter)

    df_assim_norm_v2['dt'] = [
        calendar.timegm(
            df_assim_norm_v2['dt'][i].timetuple()
        )
        for i in range(len(df_assim_norm_v2['dt']))
    ]

    df_ol['dt'] = [calendar.timegm(df_ol['dt'][i].timetuple()) for i in range(len(df_ol['dt']))]

    df_ol.sort_values(by=['dt'], inplace=True)

    df_meas['dt'] = [calendar.timegm(df_meas['dt'][i].timetuple()) + cst_cdt_con(df_meas['tz'][i]) for i in
                     range(len(df_meas['dt']))]

    df_meas['Q'] = pd.to_numeric(df_meas['Q'], errors='coerce')
    df_meas.dropna(inplace=True)
    df_meas['Q'] *= cfs2cms

    data = dict()
    df_names = ['Q', 'Q_ol', 'Q_assim_norm_v2']
    data['Q'] = df_meas
    data['Q_ol'] = df_ol
    data['Q_assim_norm_v2'] = df_assim_norm_v2
    for simulation in ['Q', 'Q_ol', 'Q_assim_norm_v2']:
        data[simulation]['dt'] = [datetime.utcfromtimestamp(data[simulation]['dt'][i]) for i in
                                  range(len(data[simulation]['dt']))]
    for simulation in ['Q', 'Q_ol', 'Q_assim_norm_v2']:
        for year in [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]:
            folder_out = path.join(out_folder, str(year), simulation)
            fn_out = path.join(folder_out, str(lid) + '.csv')
            if not path.exists(folder_out):
                os.makedirs(folder_out)
            if path.exists(fn_out):
                os.remove(fn_out)
            result = data[simulation]
            if len(result) == 0:
                continue

            subset_year = result[result['dt'].dt.year == year]
            subset_year.to_csv(fn_out, float_format='%.3f', columns=['dt', 'Q'], index=False)


# %%
# LOOP OVER THE GAGES INSIDE THE GAGE LIST: lid_table
# lid: LINK_ID OF THE MODEL
#
lid_list = lid_table['lid']
# t_0 = datetime.now()
# # Multithread computation starts
n_cpu = cpu_count()
pool = ThreadPool(n_cpu - 1)
results = pool.map(calc_metrics, lid_list)
# # Multithread computation ends
# t_1 = datetime.now()


#%%
# lid_list = lid_table['lid']
# lid = 230601
# calc_metrics(lid)
