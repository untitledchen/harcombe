import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pdb

data_mono = pd.read_csv('times_mono_450Ta3.csv', na_filter=False)
data_co = pd.read_csv('times_co_215Ta3.csv', na_filter=False)
#input('ENTER CO file_name') + '.csv'

data = pd.concat([data_mono, data_co], axis=0, ignore_index=True)

#sns.relplot(data=data, kind='line', x='gen', y='tol_time', hue='rep', style='culture', dashes=False, markers=True)

# data_co_init1010 = pd.read_csv('times_co_248Ta3init1010.csv', na_filter=False)
# data_co_init1010 = data_co_init1010.assign(culture='coinit')
#
# data_co_sstag = pd.read_csv('times_co_971Ta3sstag.csv', na_filter=False)
# data_co_sstag = data_co_sstag.assign(culture='cosstag')
#
# data = pd.concat([data_mono, data_co, data_co_init1010, data_co_sstag], axis=0, ignore_index=True)
# sns.relplot(data=data, kind='line', x='gen', y='tol_time', hue='culture', style='culture', dashes=False, markers=True)

slopemono, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data['culture'] == 'mono']['cycle'], data.loc[data['culture'] == 'mono']['tol_time'])
slopeco, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data['culture'] == 'co']['cycle'], data.loc[data['culture'] == 'co']['tol_time'])
sns.lmplot(data=data, x='cycle', y='tol_time', hue='culture')
print(slopemono)
print(slopeco)
plt.show()