import pdb
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from calc_frac import run_calc_frac

seed = random.randrange(1000)
phase0_lengths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
def generate():
    from results_29062023.frac_slope_over_phase0_duration.hcb_sim_phase0 import run

    for phase0_length in phase0_lengths:
        run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), phase0_length)
        run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), phase0_length)

def graph():
    fracs_mono = []
    fracs_co = []
    for phase0_length in phase0_lengths:
        last_cyc_mono = run_calc_frac(f'hcb_sim_mono_166_met1000_phase0{phase0_length}.csv', 1000, 5, last_cyc_only=True, file_write=False)
        fracs_mono.append(last_cyc_mono['frac'].mean())

        last_cyc_co = run_calc_frac(f'hcb_sim_co_166_met1_phase0{phase0_length}.csv', 1000, 5, last_cyc_only=True, file_write=False)
        fracs_co.append(last_cyc_co['frac'].mean())

    plt.plot(phase0_lengths, fracs_mono, label="mono")
    plt.plot(phase0_lengths, fracs_co, label="co")
    plt.legend()
    plt.show()

generate()