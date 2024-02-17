import pdb
import random

import sys
sys.path.insert(0, "C:\\Users\\untit\\harcombe")

from onewayanova import run_OWA

seed = 245 #random.randrange(1000)
def generate_new():
    from hcb_sim import run

    run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
    run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

    return

def calc_new():
    from calc_frac import run_calc_frac
    run_calc_frac(f'C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', 1000, 5)
    run_calc_frac(f'C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', 1000, 5)

    return

def analyze_OWA():
    run_OWA(f'frac_hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', "frac")

    return

def graph():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    # def rand_jitter(arr):
    #     stdev = .01 * (max(arr) - min(arr))
    #     return arr + np.random.randn(len(arr)) * stdev

    data1 = pd.read_csv(f'frac_hcb_sim_{"mono"}_{seed}_met{1000}_lac{1000}.csv', header=2)
    data1_groupby_cyc = data1.groupby(['cycle'])['frac'].agg(["mean", "std"])
    data2 = pd.read_csv(f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', header=2)
    data2_groupby_cyc = data2.groupby(['cycle'])['frac'].agg(["mean", "std"])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.scatter(data1['cycle'], data1['frac'], facecolor="#41A1F8", s=10, linewidths=1, label="Monoculture")
    ax.errorbar(data1_groupby_cyc.index.to_numpy(), data1_groupby_cyc["mean"].values, yerr=data1_groupby_cyc["std"].values, color="#41A1F8", linewidth=2, linestyle="--", ecolor="#054f94", elinewidth=1, capsize=4)

    offset = 0.15
    ax.scatter(data2['cycle'] + offset, data2['frac'], facecolor="#80D554", s=10, linewidths=1, label="Coculture")
    ax.errorbar(data2_groupby_cyc.index.to_numpy() + offset, data2_groupby_cyc["mean"].values, yerr=data2_groupby_cyc["std"].values, color="#80D554", linewidth=2, linestyle="--", ecolor="#479023", elinewidth=1, capsize=4)

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Tolerance (5-hour Survival Fraction)")
    ax.grid(linestyle="--", axis='y')
    ax.tick_params(direction='in')
    ax.set_xticks(np.arange(0, 10, step=1))  # Set label locations.
    ax.set_xticklabels(np.arange(1, 11, step=1))
    ax.set_ylim([0, 0.18])
    ax.legend(loc='lower right')

    plt.show()

    return

# generate_new()
# calc_new()
# analyze_OWA()
graph()

