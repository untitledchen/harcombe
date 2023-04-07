# "1) x-axis is the starting methionine concentration in monoculture
# (co-culture methionine is always 1), y-axis is lactose concentration in both mono- and co-culture.
# Color denotes ratio between evolutionary rate of tolerance."

# 2) x-axis is the starting methionine concentration in monoculture
# (co-culture methionine is always 1), y-axis is lactose concentration in both mono- and co-culture.
# Color denotes ratio between evolutionary rate of 5-hour survival fraction.

# 3) Keep lactose concentration 1000 (which fits my PDE model setup and all other modeling work
# on the system in the Harcombe lab). x-axis is starting methionine concentration in monoculture
# (co-culture methionine is always 1), y-axis is the E. coli : S. enterica growth rate ratios.
# Color denotes ratio between evolutionary rate of tolerance.

# 4) Keep lactose concentration 1000 (which fits my PDE model setup and all other modeling work
# on the system in the Harcombe lab). x-axis is starting methionine concentration in monoculture
# (co-culture methionine is always 1), y-axis is the E. coli : S. enterica growth rate ratios.
# Color denotes ratio between evolutionary rate of 5-hour survival fraction."

# met 1 to 1600
# lac 1 to 1600

# 8 steps

# calc tol at .1 accuracy to speed up
# you can fix tuple declaration in the analog later

import random
import time
import pandas as pd
from hcb_sim_heatmap import run
from calc_tolerance_heatmap import run_calc_tol
from calc_frac_heatmap import run_calc_frac
from evoratio import run_evoratio

# met_incs = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]
# lac_incs = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]
# #
# seed = 0
# mdkratios = [("seed", "met", "lac", "mdkratio")]
# fracratios = [("seed", "met", "lac", "fracratio")]
# for met_inc in met_incs:
#     seed += 1
#     print(seed)#
#
#     for lac_inc in lac_incs:
#         run(seed, "mono", 5, (0.0003, 0), 10, (met_inc, lac_inc, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
#         run(seed, "co", 5, (0.0003, 0.0003), 10, (1, lac_inc, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
#
#         file_mono = f"hcb_sim_mono_{seed}_met{met_inc}_lac{lac_inc}.csv"
#         file_co = f"hcb_sim_co_{seed}_met1_lac{lac_inc}.csv"
#
#         # MDK99 ratio
#         run_calc_tol(file_mono, 1000, 0.01, 0.1)
#         run_calc_tol(file_co, 1000, 0.01, 0.1)
#
#         tol_file_mono = f"tol_{file_mono}"
#         tol_file_co = f"tol_{file_co}"
#
#         mdkratio = run_evoratio(tol_file_mono, tol_file_co, "tol_time")
#
#         mdkratios.append((seed, met_inc, lac_inc, mdkratio)) # this is the ratio of the logged MDKs!!!
#
#         # frac ratio
#         run_calc_frac(file_mono, 1000, 5)
#         run_calc_frac(file_co, 1000, 5)
#
#         frac_file_mono = f"frac_{file_mono}"
#         frac_file_co = f"frac_{file_co}"
#
#         fracratio = run_evoratio(frac_file_mono, frac_file_co, "frac")
#
#         fracratios.append((seed, met_inc, lac_inc, fracratio)) # this is the ratio of the logged fracs!!!
#
# mdkratios_pd = pd.DataFrame(mdkratios[1:], columns=list(mdkratios[0]))
# mdkratios_pd.to_csv(f'mdkratios_met{met_incs}_lac{lac_incs}.csv', index=False)
#
# fracratios_pd = pd.DataFrame(fracratios[1:], columns=list(fracratios[0]))
# fracratios_pd.to_csv(f'fracratios_met{met_incs}_lac{lac_incs}.csv', index=False)

###
met_incs = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]
rs_incs = [0.25, 0.5, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0]
#
seed = 0
mdkratios = [("seed", "met", "rratio", "mdkratio")]
fracratios = [("seed", "met", "rratio", "fracratio")]

begin = time.perf_counter()#
for met_inc in met_incs:
    seed += 1
    print(seed)#

    for rs_inc in rs_incs:
        run(seed, "mono", 10, (0.0005, 0), 10, (met_inc, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), rs_inc)
        run(seed, "co", 10, (0.0005, 0.0005), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), rs_inc)

        file_mono = f"hcb_sim_mono_{seed}_met{met_inc}_rs{rs_inc}.csv"
        file_co = f"hcb_sim_co_{seed}_met1_rs{rs_inc}.csv"

        # MDK99 ratio
        run_calc_tol(file_mono, 1000, 0.01, 0.1, rs_inc)
        run_calc_tol(file_co, 1000, 0.01, 0.1, rs_inc)

        tol_file_mono = f"tol_{file_mono}"
        tol_file_co = f"tol_{file_co}"

        mdkratio = run_evoratio(tol_file_mono, tol_file_co, "tol_time")

        mdkratios.append((seed, met_inc, 1/rs_inc, mdkratio)) # this is the ratio of the logged MDKs!!!

        # frac ratio
        run_calc_frac(file_mono, 1000, 5, rs_inc)
        run_calc_frac(file_co, 1000, 5, rs_inc)

        frac_file_mono = f"frac_{file_mono}"
        frac_file_co = f"frac_{file_co}"

        fracratio = run_evoratio(frac_file_mono, frac_file_co, "frac")

        fracratios.append((seed, met_inc, 1/rs_inc, fracratio)) # this is the ratio of the logged fracs!!!

    print(time.perf_counter() - begin)  #
    begin = time.perf_counter()  #

mdkratios_pd = pd.DataFrame(mdkratios[1:], columns=list(mdkratios[0]))
mdkratios_pd.to_csv(f'mdkratios_met{met_incs}_rs{rs_incs}.csv', index=False)

fracratios_pd = pd.DataFrame(fracratios[1:], columns=list(fracratios[0]))
fracratios_pd.to_csv(f'fracratios_met{met_incs}_rs{rs_incs}.csv', index=False)
