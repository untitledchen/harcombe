import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb#

data = pd.read_csv(input(), na_filter=False).pivot('met', 'rratio', 'fracratio')
sns.heatmap(data).invert_yaxis()
plt.show()