import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb#

data = pd.read_csv(input(), na_filter=False).pivot('met', 'lac', 'mdkratio')
sns.heatmap(data)
plt.show()

# problem spots: lac 250 vs 1000 or 1500 or 2000 met