from hcb_sim import run
from calc_frac import run_calc_frac
import pandas as pd
import random
import numpy as np
from scipy.stats import linregress

mets = [1]
lacs = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 2000, 2500, 3000]

seed = 580#random.randrange(1000)

frac_rates = [("culture", "rep", "met", "lac", "frac_rate")]
for met in mets:
    print(met)#
    for lac in lacs:
        # run(seed, "mono", 10, (0.01, 0), 10, (met, lac, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
        # run(seed, "co", 10, (0.01, 0.01), 10, (met, lac, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
        #
        # # file_mono = f"C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_mono_{seed}_met{met}_lac{lac}.csv"
        # file_co = f"C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_co_{seed}_met{met}_lac{lac}.csv"
        #
        # # run_calc_frac(file_mono, 1000, 5)
        # run_calc_frac(file_co, 1000, 5)

        # frac_file_mono = f"C:\\Users\\untit\\harcombe\\data_hold\\frac_hcb_sim_mono_{seed}_met{met}_lac{lac}.csv"
        frac_file_co = f"C:\\Users\\untit\\harcombe\\data_hold\\frac_hcb_sim_co_{seed}_met{met}_lac{lac}.csv"

        # compute per rep
        # data1 = pd.read_csv(f'frac_hcb_sim_mono_685_met1000_phase0{phase}.csv', header=2, na_filter=False)
        data2 = pd.read_csv(frac_file_co, header=2, na_filter=False)

        yvar = 'frac'
        # data1[f'log10_{yvar}'] = np.log10(data1[yvar])
        data2[f'log10_{yvar}'] = np.log10(data2[yvar])

        reps = range(max(data2['rep']) + 1)
        #reps2 = range(max(data2['rep']) + 1)
        for rep in reps:   #とりあえず同じだし
            # data1_rep = data1.loc[data1['rep'] == rep]
            data2_rep = data2.loc[data2['rep'] == rep]

            # slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
            # lst.append(('Monoculture', phase, rep, slope_mono))

            slope_co, intercept, r_value, p_value, std_err = linregress(data2_rep['cycle'], data2_rep[f'log10_{yvar}'])
            frac_rates.append(('Coculture', rep, met, lac, slope_co))

frac_rates_pd = pd.DataFrame(frac_rates[1:], columns=list(frac_rates[0]))
frac_rates_pd.to_csv(f'frac_rates_met{mets}_lac{lacs}.csv', index=False)