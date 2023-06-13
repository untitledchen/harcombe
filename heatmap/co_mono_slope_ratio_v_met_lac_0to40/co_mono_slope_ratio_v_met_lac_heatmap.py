from hcb_sim import run
from calc_frac import run_calc_frac
from heatmap.co_mono_slope_ratio_v_met_lac_0to40.evoratio import run_evoratio
import pandas as pd

# mets = [0, 1, 5, 10, 15, 20]
# lacs = [0, 1, 5, 10, 15, 20]
mets = [35, 40]
lacs = [30]

seed = 339 # random.randrange(1000)

fracratios = [("seed", "met", "lac", "fracratio")]
for met in mets:
    print(met)#
    for lac in lacs:
        run(seed, "mono", 10, (0.01, 0), 10, (met, lac, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
        run(seed, "co", 10, (0.01, 0.01), 10, (met, lac, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

        file_mono = f"C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_mono_{seed}_met{met}_lac{lac}.csv"
        file_co = f"C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_co_{seed}_met{met}_lac{lac}.csv"

        run_calc_frac(file_mono, 1000, 5)
        run_calc_frac(file_co, 1000, 5)

        frac_file_mono = f"C:\\Users\\untit\\harcombe\\data_hold\\frac_hcb_sim_mono_{seed}_met{met}_lac{lac}.csv"
        frac_file_co = f"C:\\Users\\untit\\harcombe\\data_hold\\frac_hcb_sim_co_{seed}_met{met}_lac{lac}.csv"

        fracratio = run_evoratio(frac_file_mono, frac_file_co, "frac")

        fracratios.append((seed, met, lac, fracratio)) # this is the ratio of the logged fracs!!!

fracratios_pd = pd.DataFrame(fracratios[1:], columns=list(fracratios[0]))
fracratios_pd.to_csv(f'fracratios_met{mets}_lac{lacs}.csv', index=False)