import pdb

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from calc_frac import run_calc_frac

xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
fracs_mono = []
fracs_co = []
for x in xs:

    last_cyc_mono = run_calc_frac(f'hcb_sim_mono_166_met1000_phase0{x}.csv', 1000, 5, last_cyc_only=True, file_write=False)
    fracs_mono.append(last_cyc_mono['frac'].mean())
    # pdb.set_trace()

    last_cyc_co = run_calc_frac(f'hcb_sim_co_166_met1_phase0{x}.csv', 1000, 5, last_cyc_only=True, file_write=False)
    fracs_co.append(last_cyc_co['frac'].mean())

plt.plot(xs, fracs_mono, label="mono")
plt.plot(xs, fracs_co, label="co")
plt.legend()
plt.show()
