import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb#

data = pd.read_csv(input(), na_filter=False).pivot('met', 'rratio', 'fracratio')
sns.heatmap(data).invert_yaxis()
plt.show()

# problem spots: lac 250 vs 1000 or 1500 or 2000 met
#7,1750,500,75.38964371059546
#3,750,250,19.91186149668294
#2,500,750,12.818114452170775f
#, vmax=3.0