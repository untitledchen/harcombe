import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv(
    "frac_rates_met[1000]_lac[1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500].csv")
x = sns.relplot(data=data, x="lac", y="frac_rate", hue="culture", kind="line", ci="sd", err_style="bars", alpha=0.65, palette="husl")
x.set(xlabel="Initial Lactose", ylabel="Slope of Lin. Reg. Line of Log10(Survival Fraction) by Rep", title="Impact of Initial Lactose on Monoculture Survival Fraction Evolution")
plt.subplots_adjust(top=0.95) # use a lower number to make more vertical space
plt.show()