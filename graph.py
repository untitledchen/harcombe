import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb#

data = pd.read_csv(input())
x = sns.relplot(data=data, x="lac", y="frac_rate", hue="culture", kind="line", ci="sd", err_style="bars", alpha=0.65, palette="husl")
x.set(xlabel="Initial Lactose", ylabel="Slope of Lin. Reg. Line of Log10(Survival Fraction) by Rep", title="Impact of Initial Lactose on Conoculture Survival Fraction Evolution")
plt.subplots_adjust(top=0.95) # use a lower number to make more vertical space
plt.show()