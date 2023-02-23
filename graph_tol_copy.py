import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pdb

fileco = "times_init_pop1000_perc_cutoff0.01_interval0.25_costag_co_seed134_rep20_mu(0.01, 0.01)_cycles10_init_R(1, 2780, 0)_init_n(5, 5)_init_lag(1, 1)_Ta5_alpha(3, 3)_null(1.1, 1.1).csv"
filemono = "times_init_pop1000_perc_cutoff0.01_interval0.25_final_mono_seed384_rep20_mu(0.01, 0.01)_cycles10_init_R(80, 2780, 0)_init_n(5, 5)_init_lag(1, 1)_Ta5_alpha(3, 3)_null(1.1, 1.1).csv"

data_mono = pd.read_csv(filemono, na_filter=False)
data_co = pd.read_csv(fileco, na_filter=False)

data = pd.concat([data_mono, data_co], axis=0, ignore_index=True)

# sns.relplot(data=data, kind='line', x='cycle', y='tol_time', hue='rep', style='culture', dashes=False, markers=True)

# data_co_init1010 = pd.read_csv('times_co_248Ta3init1010.csv', na_filter=False)
# data_co_init1010 = data_co_init1010.assign(culture='coinit')
#
# data_co_sstag = pd.read_csv('times_co_971Ta3sstag.csv', na_filter=False)
# data_co_sstag = data_co_sstag.assign(culture='cosstag')
#
# data = pd.concat([data_mono, data_co, data_co_init1010, data_co_sstag], axis=0, ignore_index=True)
data = pd.concat([data_mono, data_co], axis=0, ignore_index=True)
sns.relplot(data=data, kind='line', x='cycle', y='tol_time', hue='culture', style='culture', dashes=False, markers=True)
#
slopemono, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data['culture'] == 'mono']['cycle'], data.loc[data['culture'] == 'mono']['tol_time'])
slopeco, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data['culture'] == 'co']['cycle'], data.loc[data['culture'] == 'co']['tol_time'])
sns.lmplot(data=data, x='cycle', y='tol_time', hue='culture')
print(slopemono)
print(slopeco)
plt.show()