import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from tolerance_odes import odes_co

import random
import copy
import math
import pdb#

#question
## should 105 mutants start lagging?
## -lag pop in gen 2 only: -= 1 in mutation probably
# no 0.1 K_M
# check mutation round_half_up
# assuming genetic composition stays exactly the same during 500 mic transfer
# more tuples

globals()['seed'] = random.randrange(1000)
random.seed(seed) #891
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

def run_phase(init_cond, lags, t_interval, phase, names):
    alpha = [2,0][phase-1]
    sol = odeint(odes_co, init_cond, t_interval, args=(alpha, lags, False))

    if globals()['genx'] == 19:
        flatlags = [element for sub_list in lags for element in sub_list]
        flatnames = [element for sub_list in names for element in sub_list]

        j = 0
        for i in range(0, int(sol[:, 3:].shape[1]), 2):
            #pdb.set_trace()
            plt.plot(t_interval, np.add(sol[:, i + 3], sol[:, i + 4]), label=f'{round_half_up(flatlags[j], 3)}, {flatnames[j]}')
            j += 1

        plt.plot(t_interval, sol[:, 0], label='met', color='r')
        plt.plot(t_interval, sol[:, 1], label='lac', color='b')
        plt.plot(t_interval, sol[:, 2], label='ace', color='y')
        plt.legend(loc='center right')
        plt.title(f'seed {seed}')
        plt.show()

    return sol

# edited 10/17: returns genotype_n_sep edited
def generate_mutants(genotype_n_sep, nE, u, mutation_function, gen):
    genotype_n_growing = [genotype_n_sep[1:(2*nE + 1):2], genotype_n_sep[(2*nE + 1)::2]] # makes a deep copy, apparently    ##

    for species, growing in enumerate(genotype_n_growing):
        genotype_freq = [i/sum(growing) for i in growing]  ## faster to use numpy or sth?

        # if round_half_up(sum(genotype_freq)) != 1: #
        #     print(genotype_freq)
        #     print("stop freq")
        #     return -1

        # Adamowicz-based
        chance = [random.uniform(0, 1) for i in range(round_half_up(sum(growing)))]
        chance_tf = [i < u[species] for i in chance]
        mutant_n = sum(chance_tf)
        print(mutant_n, end=' ')  #

        if mutant_n != 0:
            mutants = random.choices(range(len(growing)), weights=genotype_freq, k=mutant_n)
            genotype_ct = len(flask[species].genotypes)

            for j, anc_i in enumerate(mutants):
                ancestor = flask[species].genotypes[anc_i]

                genotype_n_sep[2*anc_i + (1, 2*nE + 1)[species]] -= 1 # tuple is for correct indices in genotype_n_sep, which is not separated by species
                genotype_n_sep.append(1) ##
                genotype_n_sep.append(0)

                flask[species].add_genotype(Genotype(f"{('E', 'S')[species]}{genotype_ct + j}g{gen}", n=1, lag=max([0, ancestor.lag + mutation_function()]),
                                               ancestors=ancestor.name + ' ' + ancestor.ancestors))

    return genotype_n_sep

def run_one_simulation(flask, init_R, inher_R, Ta, rep, gen, mutation_function):
    final_sub = []

    # init_cond1
    init_cond1 = [init_R[0] + inher_R[0], init_R[1] + inher_R[1], init_R[2] + inher_R[2]]
    for species in flask:
        for genotype in species.genotypes: # for each genotype, l = n and g = 0
            init_cond1.append(genotype.n)
            init_cond1.append(0)
    lags1 = [[i.lag for i in j.genotypes] for j in flask]
    t_interval1 = np.linspace(0, Ta, 1000)
    names1 = [[i.name for i in j.genotypes] for j in flask] #


    # phase 1
    sol1 = run_phase(init_cond1, lags1, t_interval1, 1, names1)

    # collect 1
    genotype_n_sep1 = list(sol1[-1, 3:])
    #print('sep1', genotype_n_sep1)#

    # append 1
    nE = len(lags1[0])#
    genotype_n_unsep1 = [genotype_n_sep1[i] + genotype_n_sep1[i + 1] for i in range(0, len(genotype_n_sep1), 2)]
    for s, species in enumerate(flask):
        for i, genotype in enumerate(species.genotypes):
            final_sub.append((rep, gen, 1, species.name, genotype.name, genotype_n_sep1[2*i + (0, 2*nE)[s]], genotype_n_sep1[2*i + (1, 2*nE + 1)[s]], genotype_n_unsep1[i + (0, nE)[s]], genotype.lag, Ta)) ## attempt to fix

    # mutation
    genotype_n_sep_mut = generate_mutants(copy.deepcopy(genotype_n_sep1), len(flask[0].genotypes), (flask[0].u, flask[1].u), mutation_function, gen)
    #print('mut', genotype_n_sep_mut)#

    # init_cond2
    init_cond2 = list(init_R)
    for i in range(0, len(genotype_n_sep_mut), 2):
        init_cond2.append(genotype_n_sep_mut[i])
        init_cond2.append(genotype_n_sep_mut[i + 1])
    #print('init cond2', init_cond2[2:]) #
    lags2 = [[i.lag for i in j.genotypes] for j in flask]
    t_interval2 = np.linspace(0, 20, 1000) # arbitrary
    names2 = [[i.name for i in j.genotypes] for j in flask]  #

    # phase 2
    sol2 = run_phase(init_cond2, lags2, t_interval2, 2, names2)

    # collect 2
    genotype_n_sep2 = list(sol2[-1, 3:])
    #print('gen sep2', genotype_n_sep2)  #
    genotype_n_unsep2 = [genotype_n_sep2[i] + genotype_n_sep2[i+1] for i in range(0, len(genotype_n_sep2), 2)]

    # register final counts of genotypes in flask
    nE = len(flask[0].genotypes)
    for i, n in enumerate(genotype_n_unsep2[:nE]):
        flask[0].genotypes[i].n = n / 2  # divide by 2

    for j, n in enumerate(genotype_n_unsep2[nE:]):
        flask[1].genotypes[j].n = n / 2  # divide by 2

    # append 2
    for s, species in enumerate(flask):
        for i, genotype in enumerate(species.genotypes):
            final_sub.append((rep, gen, 2, species.name, genotype.name, genotype_n_sep2[2*i + (0, 2*nE)[s]]/2, genotype_n_sep2[2*i + (1, 2*nE + 1)[s]]/2, genotype_n_unsep2[i + (0, nE)[s]]/2, genotype.lag, Ta)) # half the population pipetted into the next trial

    inher_R = (sol2[-1][0]/2, sol2[-1][1]/2, sol2[-1][2]/2)
    return final_sub, inher_R

### simulation parameters
reps = 1
u = 0.001 # mutation rate
gens = 20

init_R = (1, 1000, 0) # starting (M, L, A) of each new growth flask
init_n = 10 # starting population
init_lag = 1 # starting lag
Ta = 3 # length of antibiotic treatment
max_lag_change = 1.1 # max mutation-induced lag change ## orig. antibiotic_change_per_well * 1.1

### make mutation function
# null_function grabs MIC mutations from a uniform distribution
null_function = make_null_function(max_lag_change)

### run simulation
final = [('rep', 'gen', 'phase_end', 'species', 'genotype', 'nlag', 'ngrow', 'ntot', 'lag', 'Ta')]
per_gen_data = []
for rep in range(reps):

    # set up first species
    flask = Flask()
    flask.add_species(Species('Escherichia coli', u))
    flask[0].add_genotype(Genotype('E0g0', copy.deepcopy(init_n), copy.deepcopy(init_lag), ancestors = '0')) ## lag cannot be 0 to avoid divide by zero
    final.append((rep, 0, 0, 'Escherichia coli', 'E0g0', init_n, 0, init_n, init_lag, Ta))

    flask.add_species(Species('Salmonella enterica', u))
    flask[1].add_genotype(Genotype('S0g0', copy.deepcopy(init_n), copy.deepcopy(init_lag), ancestors = '0'))
    final.append((rep, 0, 0, 'Salmonella enterica', 'S0g0', init_n, 0, init_n, init_lag, Ta))

    inher_R = (0, 0, 0)
    # run simulation
    for gen in range(gens):
        globals()['genx'] = gen
        print(genx, end=':')
        final_sub, inher_R = run_one_simulation(flask, init_R, inher_R, Ta, rep, gen, null_function)
        for row in final_sub:
            final.append(row)

final_pd = pd.DataFrame(final[1:], columns = list(final[0]))
final_pd.to_csv(f'final_co_{seed}.csv', index=False)