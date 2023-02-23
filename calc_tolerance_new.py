import pandas as pd
import numpy as np
from run_phase import run_phase
from itertools import chain, repeat

import pdb#

def calc_tolerance(init_cond, interval, lags, cutoff):
    alpha = tuple([3 for i in range(len(lags))])

    iter = 0
    while sum(init_cond[3:]) > cutoff and iter < 100:
        init_cond = run_phase(alpha, init_cond, lags, interval, 1, frid=True)
        iter += 1
        init_cond = init_cond[-1, :]
    #print(sum(init_cond[3:]))#
    return iter * interval

def run(filename, init_pop, perc_cutoff, interval):
    data = pd.read_csv(filename, na_filter=False)

    seed = data['seed'][0]
    #seed = filename.split('_')[2][:3]
    culture = data['culture'][0]
    #culture = filename.split('_')[1]

    reps = max(data['rep']) + 1
    cycles = max(data['cycle']) + 1

    p2_data = data.loc[data['phase_end'] == 2]

    times = [tuple(["seed", "culture", "rep", "cycle", "tol_time"])]
    for rep in range(reps):
        curr_rep = p2_data.loc[p2_data['rep'] == rep]
        for cycle in range(cycles):
            curr_cycle = curr_rep.loc[curr_rep['cycle'] == cycle]

            n_set = (curr_cycle['ntot'] / sum(curr_cycle['ntot'])) * init_pop # for each genotype

            init_cond = [2780, 2780, 2780]
            init_cond += list(chain.from_iterable(zip(n_set, repeat(0, len(curr_cycle))))) # alternate n(lag) and 0 for init_cond

            if culture == "mono":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag'])]
            elif culture == "co":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag']),
                        tuple(curr_cycle.loc[curr_cycle['species'] == 'Salmonella enterica']['lag'])]

            time = calc_tolerance(init_cond, interval, lags, init_pop*perc_cutoff)
            times.append((seed, culture, rep, cycle, time))

    times_pd = pd.DataFrame(times[1:], columns=list(times[0]))
    times_pd.to_csv(f'times_init_pop{init_pop}_perc_cutoff{perc_cutoff}_interval{interval}_{filename}', index=False)

run("costag_co_seed134_rep20_mu(0.01, 0.01)_cycles10_init_R(1, 2780, 0)_init_n(5, 5)_init_lag(1, 1)_Ta5_alpha(3, 3)_null(1.1, 1.1).csv", 1000, 0.01, 0.25)
#run("final_mono_seed384_rep20_mu(0.01, 0.01)_cycles10_init_R(80, 2780, 0)_init_n(5, 5)_init_lag(1, 1)_Ta5_alpha(3, 3)_null(1.1, 1.1).csv", 1000, 0.01, 0.25)
