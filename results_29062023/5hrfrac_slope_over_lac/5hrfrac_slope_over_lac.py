import pdb
import random
import sys
sys.path.insert(0, "C:\\Users\\untit\\harcombe")

seed = 92 #random.randrange(1000)

lacs = [1, 5, 10, 20, 50, 100, 125, 150, 175, 200, 250, 300, 250, 400, 500, 750, 1000, 1250, 1500, 1750, 2000]
def generate_new():
    from hcb_sim import run

    for lac in lacs:
        run(seed, "mono", 10, (0.01, 0), 10, (1000, lac, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
        run(seed, "co", 10, (0.01, 0.01), 10, (1, lac, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

    return

def calc_new():
    from calc_frac import run_calc_frac
    for lac in lacs:
        run_calc_frac(f'C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_{"mono"}_{seed}_met{1000}_lac{lac}.csv', 1000, 5)
        run_calc_frac(f'C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_{"co"}_{seed}_met{1}_lac{lac}.csv', 1000, 5)

    return

def compile_slopes(yvar):
    import pandas as pd
    import numpy as np
    from scipy.stats import linregress

    frac_slopes = [("seed", "culture", "lac", "rep", "frac_slope")]
    for lac in lacs:
        data1 = pd.read_csv(f'frac_hcb_sim_{"mono"}_{seed}_met{1000}_lac{lac}.csv', header=2)
        data2 = pd.read_csv(f'frac_hcb_sim_{"co"}_{seed}_met{1}_lac{lac}.csv', header=2)

        data1[f'log10_{yvar}'] = np.log10(data1[yvar])
        data2[f'log10_{yvar}'] = np.log10(data2[yvar])

        reps = range(max(data1['rep']) + 1)
        # reps2 = range(max(data2['rep']) + 1) #とりあえず同じだし
        for rep in reps:
            data1_rep = data1.loc[data1['rep'] == rep]
            data2_rep = data2.loc[data2['rep'] == rep]

            slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
            slope_co, intercept, r_value, p_value, std_err = linregress(data2_rep['cycle'], data2_rep[f'log10_{yvar}'])

            frac_slopes.append((seed, "Monoculture", lac, rep, slope_mono))
            frac_slopes.append((seed, "Coculture", lac, rep, slope_co))

    frac_slopes_pd = pd.DataFrame(frac_slopes[1:], columns=list(frac_slopes[0]))
    frac_slopes_pd.to_csv(f'frac_slopes_{seed}_lac{lacs[0]}to{lacs[-1]}.csv', header=True, index=False, mode="w")

    return

def graph():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import copy

    data = pd.read_csv(f'frac_slopes_{seed}_lac{lacs[0]}to{lacs[-1]}_edited.csv')
    data1 = data.loc[data['culture'] == "Monoculture"]
    data2 = data.loc[data['culture'] == "Coculture"]
    data1_groupby_lac = data1.groupby(['lac'])['frac_slope'].agg(["mean", "std"])
    data2_groupby_lac = data2.groupby(['lac'])['frac_slope'].agg(["mean", "std"])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.scatter(data1['lac'], data1['frac_slope'], facecolor="#41A1F8", s=10, linewidths=1, label="Monoculture")
    offset = 50
    ax.scatter(data2['lac'] + offset, data2['frac_slope'], facecolor="#80D554", s=10, linewidths=1,
               label="Coculture")

    ax.errorbar(data1_groupby_lac.index.to_numpy(), data1_groupby_lac["mean"].values,
                yerr=data1_groupby_lac["std"].values, color="#41A1F8", linewidth=2, linestyle="--", ecolor="#054f94",
                elinewidth=1, capsize=4)

    # below0_check = data2_groupby_lac["mean"].values - data2_groupby_lac["std"].values
    # # pdb.set_trace()
    #
    # below0 = []
    # for i in range(len(below0_check)):
    #     if below0_check[i] < 0:
    #         below0.append(i)
    #
    # rem_below0 = copy.deepcopy(data2_groupby_lac["std"].values)
    # for i in below0:
    #     rem_below0[i] = 0

    ax.errorbar(data2_groupby_lac.index.to_numpy() + offset, data2_groupby_lac["mean"].values, yerr=data2_groupby_lac["std"].values, color="#80D554", linewidth=2, linestyle="--", ecolor="#479023", elinewidth=1, capsize=4)
    # ax.errorbar([data2_groupby_lac.index.to_numpy()[i] + offset for i in below0],
    #             [data2_groupby_lac["mean"].values[i] for i in below0],
    #             yerr=[data2_groupby_lac["std"].values[i] for i in below0], ecolor="k", elinewidth=1, capsize=4, )

    ax.set_xlabel("Initial Lactose (cell equivalents/mL)")
    ax.set_ylabel("Tolerance (5-hour Survival Fraction)")
    ax.grid(linestyle="--", axis='y')
    ax.tick_params(direction='in')
    ax.legend(loc='lower right')

    plt.show()

    # data = pd.read_csv(f'frac_slopes_{seed}_lac{lacs[0]}to{lacs[-1]}.csv')
    # x = sns.relplot(data=data, x="lac", y="frac_slope", hue="culture", kind="line", ci="sd", err_style="bars",
    #                 alpha=0.7, palette="husl", legend=False)
    # x.set(xlabel="Initial Lactose", ylabel="5-hour Survival Fraction Slope")
    # plt.legend(title='Culture Type', loc="lower right", labels=['Monoculture', 'Coculture'])
    # plt.show()

    return

# generate_new()
# calc_new()
# compile_slopes("frac")
graph()