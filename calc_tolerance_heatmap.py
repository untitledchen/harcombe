import pandas as pd
from hcb_sim_heatmap import run_phase
from itertools import chain, repeat

import pdb#

def calc_tolerance(init_cond, interval, lags, cutoff):
    nE = len(lags[0])  #
    alpha = tuple([3 for i in range(len(lags))])

    iter = 0
    while sum(init_cond[3:(nE*2 + 3)]) > cutoff: #and iter < 100:
        init_cond = run_phase(alpha, init_cond, lags, interval, 1, frid=False, rs=rs) #
        iter += 1
        init_cond = init_cond[-1, :]
        #print(sum(init_cond[3:(nE*2 + 3)]))

    return iter * interval

def run_calc_tol(filename, init_pop, perc_cutoff, interval, rs): ##
    globals()['rs'] = rs ##

    with open(filename, 'r') as file:
        first_line = file.readline()
    file = open(f'tol_{filename}', 'w')
    file.write(first_line)

    data = pd.read_csv(filename, header=1, na_filter=False)

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
            times.append((seed, culture, rep, cycle, time))

    times_pd = pd.DataFrame(times[1:], columns=list(times[0]))

    file.write(f'##init_pop:{init_pop}#perc_cutoff:{perc_cutoff}#interval:{interval}\n')
    times_pd.to_csv(file, index=False, mode='a')
    #times_pd.to_csv(f'times_init_pop{init_pop}_perc_cutoff{perc_cutoff}_interval{interval}_{filename}', index=False)

#run_calc_tol(input('x'), 1000, 0.01, 0.1)
