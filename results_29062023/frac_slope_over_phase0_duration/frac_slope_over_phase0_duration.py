import pdb
import random
import pandas as pd
import matplotlib.pyplot as plt
from calc_frac import run_calc_frac

seed = 314 #random.randrange(1000)
phase0_lengths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

def test():
    from results_29062023.frac_slope_over_phase0_duration.hcb_sim_phase0 import run

    phase0_length = 1
    run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), phase0_length)
    run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), phase0_length)

    return

def generate():
    from results_29062023.frac_slope_over_phase0_duration.hcb_sim_phase0 import run

    for phase0_length in phase0_lengths:
        run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), phase0_length)
        run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), phase0_length)
        print(phase0_length)

    return

def compile_fracs():
    fracs = [('culture', 'rep', 'phase0_length', 'frac')]
    for phase0_length in phase0_lengths:
        last_cyc_mono = run_calc_frac(f'hcb_sim_mono_{seed}_met1000_phase0{phase0_length}.csv', 1000, 5, last_cyc_only=True, file_write=False)
        last_cyc_co = run_calc_frac(f'hcb_sim_co_{seed}_met1_phase0{phase0_length}.csv', 1000, 5, last_cyc_only=True, file_write=False)
        reps = range(max(last_cyc_mono['rep']) + 1)

        # pdb.set_trace()
        for rep in reps:
            fracs.append(("Monoculture", rep, phase0_length, last_cyc_mono.loc[last_cyc_mono['rep'] == rep, 'frac'].item()))
            fracs.append(("Coculture", rep, phase0_length, last_cyc_co.loc[last_cyc_co['rep'] == rep, 'frac'].item()))

    frac_slopes_pd = pd.DataFrame(fracs[1:], columns=list(fracs[0]))
    frac_slopes_pd.to_csv(f'fracs_{seed}_phase0.csv', header=True, index=False, mode="w")

    return

def graph():
    import seaborn as sns

    data = pd.read_csv(f'fracs_{seed}_phase0.csv')
    x = sns.relplot(data=data, x="phase0_length", y="frac", kind="line", hue="culture", ci="sd", err_style="bars",
                    alpha=0.7, palette="husl", legend=False)
    x.set(xlabel="Phase 0 Length (hr)", ylabel="5-hour Survival Fraction")
    plt.legend(title='Culture Type', loc="lower right", labels=['Monoculture', 'Coculture'])
    plt.show()

    return

# test()
# generate()
compile_fracs()
graph()