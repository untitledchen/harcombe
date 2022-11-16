import math
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tolerance_odes import odes_mono, odes_co

from itertools import chain, repeat

import pdb

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

def odes_mono(init_cond, t_interval, a, lags):
    nE = len(lags)

    for i in range(nE):
        locals()[f'El{i}'] = init_cond[2*i]
        locals()[f'Eg{i}'] = init_cond[2*i + 1]

        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/lags[i]
        locals()[f'dEg{i}dt'] = a*locals()[f'Eg{i}'] + locals()[f'El{i}']/lags[i]

    to_return = []
    for i in range(nE):
        to_return.append(locals()[f'dEl{i}dt'])
        to_return.append(locals()[f'dEg{i}dt'])

    return to_return

def odes_co(init_cond, t_interval, a, lags):
    nE = len(lags[0])
    nS = len(lags[1])

    for i in range(nE):
        locals()[f'El{i}'] = init_cond[2*i]
        locals()[f'Eg{i}'] = init_cond[2*i + 1]

        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/lags[0][i]
        locals()[f'dEg{i}dt'] = a*locals()[f'Eg{i}'] + locals()[f'El{i}']/lags[0][i]

    for j in range(nS):
        locals()[f'Sl{j}'] = init_cond[2*j + nE]
        locals()[f'Sg{j}'] = init_cond[2*j + 1 + nE]

        locals()[f'dSl{j}dt'] = -locals()[f'Sl{j}']/lags[1][j]
        locals()[f'dSg{j}dt'] = a*locals()[f'Sg{j}'] + locals()[f'Sl{j}']/lags[1][j]

    to_return =[]
    for i in range(nE):
        to_return.append(locals()[f'dEl{i}dt'])
        to_return.append(locals()[f'dEg{i}dt'])
    for j in range(nS):
        to_return.append(locals()[f'dSl{j}dt'])
        to_return.append(locals()[f'dSg{j}dt'])

    return to_return

def calc_tolerance(culture, init_cond, lags, t_interval, stop):
    iter = 0
    sol = init_cond
    while iter < 100:
        sol = odeint([odes_mono, odes_co][culture], sol, t_interval, args=(-1, lags))[-1]

        if iter == 0:
            plt.plot(sol, t_interval)
            plt.show()

        if np.sum(sol[-1]) <= stop:
            return iter
        iter += 1
    return -1

### start interface
file_name = 'final_co_215Ta3.csv' #'ENTER FILENAME') + '.csv'
culture_type = file_name.split('_')[1]
seed = file_name.split('_')[2].split('.')[0]

data = pd.read_csv(file_name, na_filter=False)
set_pop = 1000
one_perc = set_pop * 0.01
interval = 1
n_species = [0, 1][culture_type == 'co']
reps = max(data['rep']) + 1
gens = max(data['gen']) + 1

p2_data = data.loc[data['phase_end'] == 2]

t_interval = np.linspace(0, interval, 1000)

times = [tuple(["seed", "culture", "rep", "gen", "tol_time"])]
for rep in range(reps):
    curr_rep = p2_data.loc[p2_data['rep'] == rep]
    for gen in range(gens):
        curr_gen = curr_rep.loc[curr_rep['gen'] == gen]

        n_set = (curr_gen['ntot'] / sum(curr_gen['ntot'])) * set_pop

        init_cond = list(chain.from_iterable(zip(n_set, repeat(0, len(curr_gen)))))

        if n_species == 0:
            lags = tuple(curr_gen.loc[curr_gen['species'] == 'Escherichia coli']['lag'])
            time = calc_tolerance(0, copy.deepcopy(init_cond), lags, t_interval, one_perc)
        elif n_species == 1:
            lags = [tuple(curr_gen.loc[curr_gen['species'] == 'Escherichia coli']['lag']),
                    tuple(curr_gen.loc[curr_gen['species'] == 'Salmonella enterica']['lag'])]
            time = calc_tolerance(1, copy.deepcopy(init_cond), lags, t_interval, one_perc)

        times.append((seed, ["mono", "co"][n_species], rep, gen, time))

times_pd = pd.DataFrame(times[1:], columns=list(times[0]))
times_pd.to_csv(f'weirdtimes_{["mono", "co"][n_species]}_{seed}.csv', index=False)