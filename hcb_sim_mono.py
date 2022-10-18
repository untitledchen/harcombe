import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colormap

from scipy.integrate import odeint
from tolerance_odes import odes_mono
import seaborn as sns

import random
import copy
import math
import pdb#

#question
# should i start with 0s in init_R?
# check mutation round_half_up
# resources 500 mics/1 ml, but genetic composition won't change
# more tuples

globals()['seed'] = random.randrange(1000)
random.seed(seed) #804 is odd #967 w gen 50 and ta 4
print('seed:', seed)

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors) #grabbed from stackoverflow

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

### mutation function
def make_null_function(max_lag_change):
    def null_function():
        lag_change = random.uniform(0,1) * max_lag_change
        return lag_change
    return null_function
            
### model functions
class Genotype():
    def __init__(self, name, n, lag=1, ancestors='0'):
        self.name = name
        self.n = n
        self.lag = lag
        self.ancestors = ancestors

    def __str__(self):
        return "{'name': " + str(self.name) + ", 'n': " + str(self.n) + ", 'lag': " + str(self.lag) + ", 'ancestors': " + str(self.ancestors) +"}"

class Species():
    def __init__(self, name, u):
        self.name = name
        self.u = u
        self.genotypes = []

    def __str__(self):
        return "{'genotypes': " + str([i.name for i in self.genotypes]) + ", 'u': " + str(self.u) + "}"

    def add_genotype(self, genotype):
        self.genotypes.append(genotype)
    
class Flask(list):
    def __init__(self):
        list.__init__(self)

    def __str__(self):
        return str([i.name for i in self])

    def add_species(self, species):
        self.append(species)

def run_phase(init_cond, lags, t_limit, phase, gens_info):
    alpha = [2,0][phase-1]
    t_interval = np.linspace(0, t_limit, 1000)
    sol = odeint(odes_mono, init_cond, t_interval, args=(alpha, lags, False))

    if globals()['genx'] == 19:
        j = 0
        for i in range(0, int(sol[:, 2:].shape[1]), 2):
            plt.plot(t_interval, np.add(sol[:, i + 2], sol[:, i + 3]), label= f'{round_half_up(lags[j], 3)}, gen:{gens_info[j]}') # total pop
            j += 1
        plt.plot(t_interval, sol[:, 0], label='met', color='r')
        plt.plot(t_interval, sol[:, 1], label='lac', color='b')
        plt.legend(loc='center right')
        plt.title(f'seed {seed}')
        plt.show()

    return sol

# edited 10/17: returns genotype_n_sep edited
def generate_mutants(genotype_n_sep, u, mutation_function, gen):
    genotype_n_growing = genotype_n_sep[1::2] # makes a deep copy, apparently
    genotype_freq = [i / sum(genotype_n_growing) for i in genotype_n_growing]  ## faster to use numpy or sth?

    if round_half_up(sum(genotype_freq)) != 1: #
        print(genotype_freq)
        print("stop freq")
        return -1

    # Adamowicz-based
    chance = [random.uniform(0, 1) for i in range(round_half_up(sum(genotype_n_growing)))]
    chance_tf = [i < u for i in chance]
    mutant_n = sum(chance_tf)
    print(mutant_n, end=' ')  #

    if mutant_n != 0:
        mutants = random.choices(range(len(genotype_n_growing)), weights=genotype_freq, k=mutant_n)

        genotype_ct = len(flask[0].genotypes)
        for i, anc_i in enumerate(mutants):
            ancestor = flask[0].genotypes[anc_i]

            genotype_n_sep[2*anc_i + 1] -= 1
            genotype_n_sep.append(1)
            genotype_n_sep.append(0)

            flask[0].add_genotype(Genotype(f'E{genotype_ct + i}g{gen}', n=1, lag=max([0, ancestor.lag + mutation_function()]),
                                           ancestors=ancestor.name + ' ' + ancestor.ancestors))

    return genotype_n_sep

def run_one_simulation(flask, init_R, Ta, rep, gen, mutation_function):
    final_sub = []

    # init_cond1
    init_cond1 = list(init_R)
    for genotype in flask[0].genotypes: # for each genotype, El = n and Eg = 0
        init_cond1.append(genotype.n)
        init_cond1.append(0)
    lags1 = [i.lag for i in flask[0].genotypes]
    gens_info1 = [i.name.split('g')[1] for i in flask[0].genotypes] #

    # phase 1
    sol1 = run_phase(init_cond1, lags1, Ta, 1, gens_info1)

    # collect 1
    #pdb.set_trace()
    genotype_n_sep1 = list(sol1[-1, 2:])
    print('gen sep1', genotype_n_sep1) #

    # append 1
    genotype_n_unsep1 = [genotype_n_sep1[i] + genotype_n_sep1[i + 1] for i in range(0, len(genotype_n_sep1), 2)]
    for i, genotype in enumerate(flask[0].genotypes):
        final_sub.append((rep, gen, 1, flask[0].name, genotype.name, genotype_n_unsep1[i], genotype.lag, Ta))

    # mutation
    genotype_n_sep_mut = generate_mutants(copy.deepcopy(genotype_n_sep1), flask[0].u, mutation_function, gen)
    #print('gen sup mut', genotype_n_sep_mut) #

    # init_cond2
    init_cond2 = list(init_R)
    for i in range(0, len(genotype_n_sep_mut), 2):
        init_cond2.append(genotype_n_sep_mut[i])
        init_cond2.append(genotype_n_sep_mut[i + 1])
    print('init cond2', init_cond2) #
    lags2 = [i.lag for i in flask[0].genotypes]
    gens_info2 = [i.name.split('g')[1] for i in flask[0].genotypes] #

    # phase 2
    sol2 = run_phase(init_cond2, lags2, 20, 2, gens_info2)

    # collect 2
    genotype_n_sep2 = list(sol2[-1, 2:])
    print('gen sep2', genotype_n_sep2)  #
    genotype_n_unsep2 = [genotype_n_sep2[i] + genotype_n_sep2[i+1] for i in range(0, len(genotype_n_sep2), 2)]

    # register final counts of genotypes in flask
    for i in range(len(genotype_n_unsep2)):
        flask[0].genotypes[i].n = genotype_n_unsep2[i]

    # append 2
    for genotype in flask[0].genotypes:
        final_sub.append((rep, gen, 2, flask[0].name, genotype.name, genotype.n/2, genotype.lag, Ta)) # half the population pipetted into the next trial

    #print('end r', sol2[-1, :2])
    init_R_new = init_R #
    #init_R_new = (init_R[0] + sol2[-1][0]/2, init_R[1] + sol2[-1][1]/2)
    return final_sub, init_R_new

### simulation parameters
reps = 1
u = 0.001 # mutation rate
gens = 20

init_R = (1000, 1000) # starting (M, L) of each new growth flask
init_n = 10 # starting E. coli population
init_lag = 1 # starting E. coli lag
Ta = 3 # length of anibiotic treatment
max_lag_change = 1.1 # max mutation-induced lag change ## orig. antibiotic_change_per_well * 1.1

### make mutation function
# null_function grabs MIC mutations from a uniform distribution
null_function = make_null_function(max_lag_change)

### run simulation
final = [('rep', 'gen', 'phase_end', 'species', 'genotype', 'n', 'lag', 'Ta')]
per_gen_data = []
for rep in range(reps):

    # set up first species
    flask = Flask()
    flask.add_species(Species('Escherichia coli', u))
    flask[0].add_genotype(Genotype('E0g0', init_n, init_lag, ancestors = '0')) ## lag cannot be 0 to avoid divide by zero

    final.append((rep, 0, 0, 'Escherichia coli', 'E0g0', init_n, init_lag, '0'))

    # run simulation
    for gen in range(gens):
        globals()['genx'] = gen
        final_sub, init_R = run_one_simulation(flask, init_R, Ta, rep, gen, null_function)
        for row in final_sub:
            final.append(row)

#pdb.set_trace()
final_pd = pd.DataFrame(final[1:], columns = list(final[0]))
final_pd.to_csv('final_mono.csv', index=False)

### plot

'''
all_data = pd.read_csv('all_data.csv', na_filter=False)

tol_data = all_data[['well', 'rep', 'gens', 'u', 'n_species', 'season', 'mutant_function']].loc[all_data['species'] == 0].loc[all_data['alive']==True].copy(deep=True)
tol_data = tol_data.groupby(['rep', 'gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).max()
tol_data = tol_data.assign(tolerance = tol_data['well'])
tol_stats = tol_data.groupby(['gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).tolerance.agg(['mean', 'std', 'count'])
tol_data = tol_data.groupby(['gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).count()[['gens', 'u', 'n_species', 'season', 'mutant_function']]
tol_data = tol_data.assign(tolerance_sd = tol_stats['std'].reset_index(drop=True).copy(deep=True), n = tol_stats['count'].reset_index(drop=True).copy(deep=True), tolerance = tol_stats['mean'].reset_index(drop=True).copy(deep=True))

error1 = tol_data['tolerance_sd'].loc[tol_data['n_species']==1] / [math.sqrt(i-1) for i in tol_data['n'].loc[tol_data['n_species']==1]]
error2 = tol_data['tolerance_sd'].loc[tol_data['n_species']==2] / [math.sqrt(i-1) for i in tol_data['n'].loc[tol_data['n_species']==2]]
error3 = tol_data['tolerance_sd'].loc[tol_data['n_species']==3] / [math.sqrt(i-1) for i in tol_data['n'].loc[tol_data['n_species']==3]]

plt.errorbar(tol_data['season'].loc[tol_data['n_species']==1], tol_data['tolerance'].loc[tol_data['n_species']==1], error1, label = '1', color = 'tab:blue')#
plt.errorbar(tol_data['season'].loc[tol_data['n_species']==2], tol_data['tolerance'].loc[tol_data['n_species']==2], error2, label = '2', color = 'tab:orange')#
plt.errorbar(tol_data['season'].loc[tol_data['n_species']==3], tol_data['tolerance'].loc[tol_data['n_species']==3], error3, label = '3', color = 'tab:green')#

plt.xlabel('transfer')
plt.ylabel('tolerance (arbitrary)')
plt.title(f'seed: {seed}')
plt.legend()
plt.show()
'''