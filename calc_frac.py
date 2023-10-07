import pandas as pd
from hcb_sim import run_phase
from itertools import repeat


def calc_frac(init_cond, duration, lags, init_pop_E):
    nE = len(lags[0])
    alpha = tuple([3 for i in range(len(lags))])
    sol = run_phase(alpha, init_cond, lags, duration, 1, frid=False)
    cells = sol[-1, 3:(nE*2 + 3)]
    frac = sum(cells) / init_pop_E
    return frac

def run_calc_frac(filepath, init_pop, duration, last_cyc_only=False, file_write=True):
    data = pd.read_csv(filepath, header=1, na_filter=False)

    culture = data['culture'][0]

    reps = max(data['rep']) + 1
    cycles = max(data['cycle']) + 1

    data['ntot'] = data['ngrow'] + data['nlag']
    p2_data = data.loc[data['phase_end'] == 2]

    fracs = [tuple(["culture", "rep", "cycle", "frac", "duration"])]

    set_cycle = 0
    if last_cyc_only:
        set_cycle = cycles - 1
    for rep in range(reps):
        curr_rep = p2_data.loc[p2_data['rep'] == rep]
        cycle = set_cycle
        while cycle < cycles:
            curr_cycle = curr_rep.loc[curr_rep['cycle'] == cycle]
            nE = len(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli'].index)
            nS = len(curr_cycle.index) - nE

            n_E = curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['ntot'] / sum(curr_cycle['ntot'])
            n_set_E = n_E * init_pop

            E_frac = sum(n_E)
            init_pop_E = init_pop * E_frac

            init_cond = [2780, 2780, 2780]
            init_cond.extend(list(repeat(0, nE)))
            init_cond.extend(list(n_set_E))

            if nS == 0:
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag'])]
            elif nS > 0:
                lags = [tuple(curr_cycle.loc[curr_cycle['species'] == 'Escherichia coli']['lag']),
                        tuple(curr_cycle.loc[curr_cycle['species'] == 'Salmonella enterica']['lag'])]
                n_set_S = curr_cycle.loc[curr_cycle['species'] == 'Salmonella enterica']['ntot'] / sum(curr_cycle['ntot']) * init_pop
                init_cond.extend(list(repeat(0, nS)))
                init_cond.extend(list(n_set_S))


            frac = calc_frac(init_cond, duration, lags, init_pop_E)
            fracs.append((culture, rep, cycle, frac, duration))

            cycle += 1

    fracs_pd = pd.DataFrame(fracs[1:], columns=list(fracs[0]))

    if file_write:
        with open(filepath, 'r') as f:
            first_line = f.readline()

        #
        filename = filepath.split("\\")[-1]

        file = open(f'frac_{filename}', 'w') #
        file.write(first_line)
        file.write(f'##init_pop:{init_pop}#duration:{duration}\n')
        fracs_pd.to_csv(file, index=False, mode='a')
        return

    return fracs_pd

#run_calc_frac(input("enter"), 1000, 5)

#x = run_calc_frac('phase0/hcb_sim_co_166_met1_phase015.csv', 1000, 5, last_cyc_only=True)
# pdb.set_trace()
