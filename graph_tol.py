import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb

data_mono = pd.read_csv('times_mono_539.csv', na_filter=False)
data_co = pd.read_csv('times_co_635.csv', na_filter=False)

data = pd.concat([data_mono, data_co], axis=0, ignore_index=True)

sns.relplot(data=data, kind='line', x='gen', y='tol_time', hue='culture', style='culture', dashes=False, markers=True)
plt.show()