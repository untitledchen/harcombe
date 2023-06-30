import pandas as pd
from hcb_sim import run_phase
from itertools import chain, repeat

import pdb#

def calc_tolerance(init_cond_now, interval, lags, cutoff):
    nE = len(lags[0])  #
    alpha = tuple([3 for i in range(len(lags))])

    iter = 0
    while True:
        init_cond_next = run_phase(alpha, init_cond_now, lags, interval*10, 1, frid=False)

        if sum(init_cond_next[-1][3:(nE*2 + 3)]) <= cutoff:
            break
        iter += 10
        init_cond_now = init_cond_next[-1, :]

    while sum(init_cond_now[3:(nE*2 + 3)]) > cutoff: #and iter < 100:
        init_cond_now = run_phase(alpha, init_cond_now, lags, interval, 1, frid=False)
        iter += 1
        init_cond_now = init_cond_now[-1, :]
        #print(sum(init_cond[3:(nE*2 + 3)]))

    return iter * interval

def run_calc_tol(filepath, init_pop, perc_cutoff, interval): ##
    with open(filepath, 'r') as f:
        first_line = f.readline()

    filename = filepath.split("\\")[-1]
    file = open(f'tol_{filename}', 'w')
    file.write(first_line)

    data = pd.read_csv(filepath, header=1, na_filter=False)

    culture = data['culture'][0]
    #culture = filename.split('_')[1]

    reps = max(data['rep']) + 1
    cycles = max(data['cycle']) + 1

    p2_data = data.loc[data['phase_end'] == 2]

    times = [tuple(["culture", "rep", "cycle", "tol_time"])]
    for rep in range(reps):
        curr_rep = p2_data.loc[p2_data['rep'] == rep]
        for cycle in range(cycles):
            curr_cycle = curr_rep.loc[curr_rep['cycle'] == cycle]

            curr_cycle = curr_cycle.assign(ntot=curr_cycle['ngrow'] + curr_cycle['nlag'])  ##
            n_set = (curr_cycle['ntot'] / sum(curr_cycle['ntot'])) * init_pop # for each genotype

            E_frac = sum(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['ntot'])  / sum(curr_cycle['ntot'])
            init_pop_E = init_pop * E_frac

            init_cond = [2780, 2780, 2780]
            init_cond += list(chain.from_iterable(zip(n_set, repeat(0, len(curr_cycle))))) # alternate n(lag) and 0 for init_cond

            if culture == "mono":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag'])]
            elif culture == "co":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag']),
                        tuple(curr_cycle.loc[curr_cycle['species'] == 'Salmonella enterica']['lag'])]

            time = calc_tolerance(init_cond, interval, lags, init_pop_E*perc_cutoff)
            times.append((culture, rep, cycle, time))

    times_pd = pd.DataFrame(times[1:], columns=list(times[0]))

    file.write(f'##init_pop:{init_pop}#perc_cutoff:{perc_cutoff}#interval:{interval}\n')
    times_pd.to_csv(file, index=False, mode='a')

# run_calc_tol(input('input'), 1000, 0.01, 0.1, 0.5)