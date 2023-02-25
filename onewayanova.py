import pdb

import pandas as pd
from scipy.stats import f_oneway

data1 = pd.read_csv(input("INPUT FILENAME 1 "), na_filter=False)
data2 = pd.read_csv(input("INPUT FILENAME 2 "), na_filter=False)

tol_times1 = list(data1['tol_time'])
tol_times2 = list(data2['tol_time'])

print(f_oneway(tol_times1, tol_times2))
###
reps1 = max(data1['rep']) + 1
cycles1 = max(data1['cycle']) + 1
avg_tol_times_by_cycle1 = []
for curr_cyc in range(cycles1):
    avg_tol_times_by_cycle1.append(sum(data1.loc[data1['cycle'] == curr_cyc]['tol_time']) / reps1)

reps2 = max(data2['rep']) + 1
cycles2 = max(data2['cycle']) + 1
avg_tol_times_by_cycle2 = []
for curr_cyc2 in range(cycles2):
    avg_tol_times_by_cycle2.append(sum(data2.loc[data2['cycle'] == curr_cyc2]['tol_time']) / reps2)

print(avg_tol_times_by_cycle1)
print(avg_tol_times_by_cycle2)
print(f_oneway(avg_tol_times_by_cycle1, avg_tol_times_by_cycle2))