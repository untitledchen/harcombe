import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb

data = pd.read_csv(input("INPUT FILENAME "), na_filter=False)
data['phase_cyc'] = data['cycle'] + data['phase_end']
data_r = pd.melt(data, id_vars=['phase_cyc'], value_vars=['M', 'L', 'A'])

sns.relplot(data=data_r, kind='line', x='phase_cyc', y='value', hue='variable', style='variable', dashes=False, markers=True)
plt.show()
pdb.set_trace()
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