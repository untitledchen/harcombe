import pdb

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, linregress

data1 = pd.read_csv(input("INPUT MONO "), na_filter=False)
data2 = pd.read_csv(input("INPUT CO "), na_filter=False)

# tol_times1 = list(data1['tol_time'])
# tol_times2 = list(data2['tol_time'])
###
reps1 = max(data1['rep']) + 1
cycles1 = max(data1['cycle']) + 1
data1['log10_tol_time'] = np.log10(data1['tol_time'])
best_fit_by_rep1 = []
for rep in range(reps1):
    slope, intercept, r_value, p_value, std_err = linregress(data1.loc[data1['rep'] == rep]['cycle'], data1.loc[data1['rep'] == rep]['log10_tol_time'])
    best_fit_by_rep1.append(slope)

reps2 = max(data2['rep']) + 1
cycles2 = max(data2['cycle']) + 1
data2['log10_tol_time'] = np.log10(data2['tol_time'])
best_fit_by_rep2 = []
for rep in range(reps2):
    slope, intercept, r_value, p_value, std_err = linregress(data2.loc[data1['rep'] == rep]['cycle'], data2.loc[data1['rep'] == rep]['log10_tol_time'])
    best_fit_by_rep2.append(slope)

slopemono, intercept, r_value, p_value, std_err = linregress(data1['cycle'], data1['log10_tol_time'])
slopeco, intercept, r_value, p_value, std_err = linregress(data2['cycle'], data2['log10_tol_time'])
print("slope mono: ", slopemono)
print("slope co: ", slopeco)

print(f_oneway(best_fit_by_rep1, best_fit_by_rep2))

# avg_tol_times_by_cycle1 = []
# for curr_cyc in range(cycles1):
#     avg_tol_times_by_cycle1.append(sum(data1.loc[data1['cycle'] == curr_cyc]['tol_time']) / reps1)
#
# reps2 = max(data2['rep']) + 1
# cycles2 = max(data2['cycle']) + 1
# avg_tol_times_by_cycle2 = []
# for curr_cyc2 in range(cycles2):
#     avg_tol_times_by_cycle2.append(sum(data2.loc[data2['cycle'] == curr_cyc2]['tol_time']) / reps2)
#
# print(avg_tol_times_by_cycle1)
# print(avg_tol_times_by_cycle2)
# print(f_oneway(avg_tol_times_by_cycle1, avg_tol_times_by_cycle2))