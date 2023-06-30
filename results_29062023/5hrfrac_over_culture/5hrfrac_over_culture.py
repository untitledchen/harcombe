import random

import sys
sys.path.insert(0, "C:\\Users\\untit\\harcombe")

from onewayanova import run_OWA

seed = 791 #random.randrange(1000)
def generate_new():
    from hcb_sim import run

    run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
    run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 10), 5, (3, 3), 42, "null", (1.1, 1.1))

def calc_new():
    from calc_frac import run_calc_frac
    run_calc_frac(f'C:\\Users\\untit\\harcombe\\results_29062023\\MDK99_over_culture\\hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', 1000, 5)
    run_calc_frac(f'C:\\Users\\untit\\harcombe\\results_29062023\\MDK99_over_culture\\hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', 1000, 5)

def analyze_OWA():
    run_OWA(f'frac_hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', "frac")

def graph():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    data1 = pd.read_csv(f'frac_hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', header=2)
    data2 = pd.read_csv(f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', header=2)
    data = pd.concat((data1, data2), ignore_index=True)
    x = sns.relplot(data=data, x="cycle", y="frac", hue="culture", kind="line", ci="sd", err_style="bars",
                    alpha=0.7, palette="husl", legend=False)
    x.set(xlabel="Cycle", ylabel="5-hour Survival Fraction")
    plt.legend(title='Culture Type', loc="lower right", labels=['Monoculture', 'Coculture'])
    plt.show()

# generate_new()
# calc_new()
# analyze_OWA()
graph()

