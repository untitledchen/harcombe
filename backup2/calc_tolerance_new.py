import pandas as pd
from backup2.run_phase import run_phase
from itertools import chain, repeat


def calc_tolerance(init_cond, interval, lags, cutoff):
    alpha = tuple([3 for i in range(len(lags))])

    iter = 0
    while sum(init_cond[3:]) > cutoff and iter < 100:
        init_cond = run_phase(alpha, init_cond, lags, interval, 1, frid=False, inc=100) #
        iter += 1

        # if iter == 1: #
        #     flatlags = [element for sub_list in lags for element in sub_list]
        #     t_interval = np.linspace(0, interval, 100)
        #     j = 0
        #     for i in range(0, int(init_cond[:, 3:].shape[1]), 2):
        #         # pdb.set_trace()
        #         plt.plot(t_interval, np.add(init_cond[:, i + 3], init_cond[:, i + 4]),
        #                  label=f'{round(flatlags[j], 3)}')
        #         j += 1
        #     plt.legend(loc='center right')
        #     plt.show()
        #     return
        init_cond = init_cond[-1, :]
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

run(input('INPUT FILENAME '), 1000, 0.01, .1)
