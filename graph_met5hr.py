from hcb_sim import run
from calc_frac import run_calc_frac

import seaborn as sns
import matplotlib.pyplot as plt
import random
import pandas as pd
import numpy as np
from scipy.stats import linregress

# seed = random.randrange(1000)
#
# lst = [('culture', 'met', 'rep', 'frac_slope')]
#
# mets = list(range(0, 10000, 500))
# for met in mets:
#     print(met)  #
#
#     run(seed, "mono", 10, (0.01, 0), 10, (met, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
#     #run(seed, "co", 5, (0.0003, 0.0003), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
#
#     run_calc_frac(f'hcb_sim_mono_{seed}_met{met}.csv', 1000, 5)
#     #run_calc_frac(f'hcb_sim_co_{seed}_met1.csv', 1000, 5)
#
#     # compute per rep
#     data1 = pd.read_csv(f'frac_hcb_sim_mono_{seed}_met{met}.csv', header=2, na_filter=False)
#     #data2 = pd.read_csv(f'frac_hcb_sim_co_{seed}_met1.csv', header=2, na_filter=False)
#
#     yvar = 'frac'
#     data1[f'log10_{yvar}'] = np.log10(data1[yvar])
#     #data2[f'log10_{yvar}'] = np.log10(data2[yvar])
#
#     reps = range(max(data1['rep']) + 1)
#     for rep in reps:
#         data1_rep = data1.loc[data1['rep'] == rep]
#         #data2_rep = data2.loc[data2['rep'] == rep]
#
#         slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
#         lst.append(('mono', met, rep, slope_mono))
#
#         # slope_co, intercept, r_value, p_value, std_err = linregress(data2_rep['cycle'], data2_rep[f'log10_{yvar}'])
#         # lst.append(('co', 1, rep, slope_co))
#
# final_pd = pd.DataFrame(lst[1:], columns=list(lst[0]))
# final_pd.to_csv(f'frac_slopes.csv', index=False)
#
# data = pd.read_csv('frac_slopes.csv')
#
# sns.relplot(data=data, x="met", y="frac_slope", kind="line", ci="sd", err_style="bars")
# plt.show()


data1 = pd.read_csv(input(), header=2, na_filter=False)
yvar = 'frac'
data1[f'log10_{yvar}'] = np.log10(data1[yvar])
slope_mono, intercept, r_value, p_value, std_err = linregress(data1['cycle'], data1[f'log10_{yvar}'])
print(slope_mono)

# frac_hcb_sim_co_561_met1.csv
# 0.11739308753272672