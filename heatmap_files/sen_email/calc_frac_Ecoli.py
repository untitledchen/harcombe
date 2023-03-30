import pandas as pd
from backup.backup2.run_phase import run_phase
from itertools import chain, repeat


def calc_frac(init_cond, duration, lags, init_pop_E):
    nE = len(lags[0])#
    alpha = tuple([3 for i in range(len(lags))])
    sol = run_phase(alpha, init_cond, lags, duration, 1, frid=False)
    cells = sol[-1, 3:(nE*2 + 3)]#
    frac = sum(cells) / init_pop_E #
    return frac

def run_calc_frac(filename, init_pop, duration):
    data = pd.read_csv(filename, na_filter=False)

    seed = data['seed'][0]
    culture = data['culture'][0]

    reps = max(data['rep']) + 1
    cycles = max(data['cycle']) + 1

    p2_data = data.loc[data['phase_end'] == 2]

    fracs = [tuple(["seed", "culture", "rep", "cycle", "frac", "duration"])]
    for rep in range(reps):
        curr_rep = p2_data.loc[p2_data['rep'] == rep]
        for cycle in range(cycles):
            curr_cycle = curr_rep.loc[curr_rep['cycle'] == cycle]

            n_set = (curr_cycle['ntot'] / sum(curr_cycle['ntot'])) * init_pop # for each genotype
            init_pop_E = sum(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['ntot']) #

            init_cond = [2780, 2780, 2780]
            init_cond += list(chain.from_iterable(zip(n_set, repeat(0, len(curr_cycle))))) # alternate n(lag) and 0 for init_cond

            if culture == "mono":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag'])]
            elif culture == "co":
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag']),
                        tuple(curr_cycle.loc[curr_cycle['species'] == 'Salmonella enterica']['lag'])]

            frac = calc_frac(init_cond, duration, lags, init_pop_E)
            fracs.append((seed, culture, rep, cycle, frac, duration))

    fracs_pd = pd.DataFrame(fracs[1:], columns=list(fracs[0]))
    fracs_pd.to_csv(f'fracs-Ecoli_init_pop{init_pop}_duration{duration}_{filename}', index=False)

run_calc_frac(input('INPUT FILENAME '), 1000, 5)
