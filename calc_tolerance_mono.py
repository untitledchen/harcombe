import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tolerance_odes import odes_mono

from itertools import chain, repeat

import pdb

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

def calc_tolerance_co(init_cond, lags, t_interval, names):
    nE = len(lags)
    alpha = 2

    sol = odeint(odes_mono, init_cond, t_interval, args=(2, lags, False))

    # half-saturation constants
    K_M = 1
    K_L = 1

    # E growth constants
    rE = 1
    kE = 5e-9

    # calc derivative - create derivatives
    derivs = []
    for t in range(t_interval.size):
        curr_sol = sol[t, :]
        curr_M = curr_sol[0]
        curr_L = curr_sol[1]

        curr_deriv = 0
        for i in range(nE):
            El = curr_sol[2 * i + 2]
            Eg = curr_sol[2 * i + 3]

            curr_deriv += -El / lags[i]
            curr_deriv += (1 - alpha) * rE * Eg * (curr_M/(curr_M + K_M)) * (curr_L/(curr_L + K_L)) - kE * Eg + El / lags[i]

        derivs.append(curr_deriv)

    # y = mx + b
    min_slope = min(derivs)
    min_ind = derivs.index(min_slope)
    min_t = t_interval[min_ind]
    y = sum(sol[min_ind, 2:])

    # plt.plot(t_interval, derivs)
    # plt.show()

    b = y - min_slope * min_t

    # plt.plot(t_interval, [np.log10(min_slope*t + b) for t in t_interval])
    # plt.show()

    # 10 = mx + b
    x = (10 - b) / min_slope
    return x

### start interface
seed = input('ENTER SEED #')
file_name = 'final_mono_' + seed + '.csv'
data = pd.read_csv(file_name, na_filter=False)
set_pop = 1000
one_perc = set_pop * 0.01
interval = 10
gens = max(data['gen']) + 1

p2_data = data.loc[data['phase_end'] == 2]

t_interval = np.linspace(0, interval, 1000)

times = [tuple(["seed", "culture", "gen", "tol_time"])]
for gen in range(gens):
    curr_gen = p2_data.loc[p2_data['gen'] == gen]

    n_set = (curr_gen['ntot'] / sum(curr_gen['ntot'])) * set_pop

    init_cond = [1000, 1000]  # ML but what should it be at?
    init_cond += list(chain.from_iterable(zip(n_set, repeat(0, len(curr_gen)))))
    lags = tuple(curr_gen.loc[curr_gen['species'] == 'Escherichia coli']['lag'])
    names = tuple(curr_gen.loc[curr_gen['species'] == 'Escherichia coli']['genotype'])

    time = calc_tolerance_co(init_cond, lags, t_interval, names)

    times.append((seed, "mono", gen, time))

times_pd = pd.DataFrame(times[1:], columns=list(times[0]))
times_pd.to_csv(f'times_mono_{seed}.csv', index=False)