import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tolerance_odes import odes_co

from itertools import chain, repeat

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

def calc_tolerance(init_cond, lags, t_interval, names):
    iter = 0
    while iter < 100:
        #sol = run_phase(init_cond, lags, t_interval, 1, names)
        sol = odeint(odes_co, init_cond, t_interval, args=(2, lags, True))
        genotype_n_sep = list(sol[-1, 3:])
        if sum(genotype_n_sep) <= one_perc:
            return interval*iter + interval
        init_cond = sol[-1, :]

        if iter%2 == 0:
            flatlags = [element for sub_list in lags for element in sub_list]
            flatnames = [element for sub_list in names for element in sub_list]

            j = 0
            for i in range(0, int(sol[:, 3:].shape[1]), 2):
                # pdb.set_trace()
                plt.plot(t_interval, np.add(sol[:, i + 3], sol[:, i + 4]),
                         label=f'{round_half_up(flatlags[j], 3)}, {flatnames[j]}')
                j += 1

            plt.plot(t_interval, sol[:, 0], label='met', color='r')
            plt.plot(t_interval, sol[:, 1], label='lac', color='b')
            plt.plot(t_interval, sol[:, 2], label='ace', color='y')
            plt.legend(loc='center right')
            #plt.title(f'seed {seed}')
            plt.show()

        iter += 1
    return sum(genotype_n_sep)

file_name = 'final_co_' + input('ENTER SEED #') + '.csv'
data = pd.read_csv(file_name, na_filter=False)
set_pop = 1000
one_perc = set_pop * 0.01
interval = 1

maxgen = max(data['gen'])
end_data = data.loc[data['gen'] == maxgen].loc[data['phase_end'] == 2]
n_set = (end_data['ntot'] / sum(end_data['ntot'])) * set_pop

init_cond = [1000, 1000, 1000] # MLA
init_cond += list(chain.from_iterable(zip(n_set, repeat(0, len(end_data)))))

lags = [tuple(end_data.loc[end_data['species'] == 'Escherichia coli']['lag']), tuple(end_data.loc[end_data['species'] == 'Salmonella enterica']['lag'])]
names = [tuple(end_data.loc[end_data['species'] == 'Escherichia coli']['genotype']), tuple(end_data.loc[end_data['species'] == 'Salmonella enterica']['genotype'])]

t_interval = np.linspace(0, interval, 1000)

print(calc_tolerance(init_cond, lags, t_interval, names))