import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb#

data = pd.read_csv(input(), na_filter=False).pivot('met', 'lac', 'fracratio')
sns.heatmap(data, vmax=10).invert_yaxis()
plt.show()