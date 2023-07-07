import random

import sys
sys.path.insert(0, "C:\\Users\\untit\\harcombe")

from onewayanova import run_OWA

seed = random.randrange(1000)

met = [1, 5, 10, 20, 50, 100, 125, 150, 175, 200, 250, 300, 250, 400, 500, 750, 1000, 1250, 1500, 1750, 2000]

def generate_new():
    from hcb_sim import run

    run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 10), 5, (3, 3), 42, "null", (1.1, 1.1))
    for met in mets:
        run(seed, "mono", 10, (0.01, 0), 10, (met, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))

    return

def calc_new():
    from calc_frac import run_calc_frac

    run_calc_frac(f'C:\\Users\\untit\\harcombe\\results_29062023\\5hrfrac_slope_over_met\\hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', 1000, 5)

    for met in mets:
        run_calc_frac(f'C:\\Users\\untit\\harcombe\\results_29062023\\5hrfrac_slope_over_met\\hcb_sim_{"mono"}_{seed}_met{met}_lac{1000}.csv', 1000, 5)

    return

def compile_slopes(yvar):
    import pandas as pd
    import numpy as np
    from scipy.stats import linregress

    # coculture
    # data2 = pd.read_csv(f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{1000}.csv', header=2)
    # data2[f'log10_{yvar}'] = np.log10(data2[yvar])
    # reps2 = range(max(data2['rep']) + 1)
    # for rep in reps2:
    #     data2_rep = data2.loc[data2['rep'] == rep]
    #     slope_co, intercept, r_value, p_value, std_err = linregress(data2_rep['cycle'], data2_rep[f'log10_{yvar}'])
    #     frac_slopes.append((seed, "Coculture", met, rep, slope_co))

    frac_slopes = [("seed", "culture", "met", "rep", "frac_slope")]
    for met in mets:
        data1 = pd.read_csv(f'frac_hcb_sim_{"mono"}_{seed}_met{met}_lac{1000}.csv', header=2)

        data1[f'log10_{yvar}'] = np.log10(data1[yvar])

        reps = range(max(data1['rep']) + 1)
        for rep in reps:
            data1_rep = data1.loc[data1['rep'] == rep]

            slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])

            frac_slopes.append((seed, "Monoculture", met, rep, slope_mono))

    frac_slopes_pd = pd.DataFrame(frac_slopes[1:], columns=list(frac_slopes[0]))
    frac_slopes_pd.to_csv(f'frac_slopes_{seed}_met{mets[0]}to{mets[-1]}.csv', header=True, index=False, mode="w")

    return

def graph():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots() ##

    data = pd.read_csv(f'frac_slopes_{seed}_met{mets[0]}to{mets[-1]}.csv')
    x = sns.relplot(data=data, x="met", y="frac_slope", kind="line", ci="sd", err_style="bars",
                    alpha=0.7, palette="husl", legend=False)
    x.set(xlabel="Initial Methionine", ylabel="5-hour Survival Fraction Slope")
    plt.legend(title='Culture Type', loc="lower right", labels=['Monoculture', 'Coculture'])
    plt.show()

    return

# generate_new()
# calc_new()
# compile_slopes("frac")
# graph()