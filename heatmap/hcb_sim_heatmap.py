import pdb
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tolerance_odes_heatmap import odes

import random
import copy
import math

# make seed visible for re-generating purposes
# globals()['seed'] = random.randrange(1000)
# random.seed(seed)
# print('seed:', seed)

##
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors)

# rounding function: rounds .5 up to next int
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

# null_function grabs MIC mutations from a uniform distribution
def make_null_function(max_lag_change):
    def null_function(s):
        lag_change = random.uniform(0,1) * max_lag_change[s]
        return lag_change
    return null_function
            
# objects
class Genotype():
    def __init__(self, name, n, lag=1, ancestors='0'):
        self.name = name
        self.n = n
        self.lag = lag
        self.ancestors = ancestors

    def __str__(self):
        return "{'name': " + str(self.name) + ", 'n': " + str(self.n) + ", 'lag': " + str(self.lag) + ", 'ancestors': " + str(self.ancestors) +"}"

class Species():
    def __init__(self, name, mu):
        self.name = name
        self.mu = mu
        self.genotypes = []

    def __str__(self):
        return "{'genotypes': " + str([i.name for i in self.genotypes]) + ", 'mu': " + str(self.mu) + "}"

    def add_genotype(self, genotype):
        self.genotypes.append(genotype)
    
class Flask(list):
    def __init__(self):
        list.__init__(self)

    def __str__(self):
        return str([i.name for i in self])

    def add_species(self, species):
        self.append(species)

def run_phase(alpha, init_cond, lags, t, phase, inc=1000, frid=False, rs=None): #
    alpha_this = tuple([[a,0][phase-1] for a in alpha])
    t_interval = np.linspace(0, t, inc)
    sol = odeint(odes, init_cond, t_interval, args=(alpha_this, lags, frid, rs)) ##

    return sol

def generate_mutants(flask, genotype_n_sep, nE, mu, mutation_function, cycle):
    genotype_n_growing = [genotype_n_sep[1:(2*nE + 1):2], genotype_n_sep[(2*nE + 1)::2]]

    for species, strain in enumerate(flask): # [genotype_n_growing[0]] to exclude S
        growing = genotype_n_growing[species]
        u = mu[species]
        genotype_freq = [i/sum(growing) for i in growing]

        # Adamowicz-based
        chance_tf = [random.uniform(0, 1) < u for i in range(round_half_up(sum(growing)))]
        mutant_n = sum(chance_tf)

        if mutant_n != 0:
            mutants = random.choices(range(len(growing)), weights=genotype_freq, k=mutant_n)
            genotype_ct = len(flask[species].genotypes)

            for j, anc_i in enumerate(mutants):
                ancestor = flask[species].genotypes[anc_i]

                genotype_n_sep[2*anc_i + (1, 2*nE + 1)[species]] -= 1 # tuple is for correct indices in genotype_n_sep, which is not separated by species
                genotype_n_sep.append(1) ##
                genotype_n_sep.append(0)

                flask[species].add_genotype(Genotype(f"{('E', 'S')[species]}{genotype_ct + j}c{cycle}", n=1, lag=max([0, ancestor.lag + mutation_function(species)]),
                                               ancestors=ancestor.name + ' ' + ancestor.ancestors))

    return genotype_n_sep

def run_one_simulation(seed, culture, flask, init_R, inher_R, Ta, alpha, t_grow, rep, cycle, mutation_function):
    final_sub = []

    # init_cond1
    init_cond = [init_R[0] + inher_R[0], init_R[1] + inher_R[1], init_R[2] + inher_R[2]]
    for species in flask:
        for genotype in species.genotypes:
            init_cond.append(genotype.n)
            init_cond.append(0)
    lags1 = [[i.lag for i in j.genotypes] for j in flask]

    # phase 1
    sol1 = run_phase(alpha, init_cond, lags1, Ta, 1, rs=rs) ##

    # collect 1
    genotype_n_sep1 = list(sol1[-1, 3:])

    # append 1
    nE = len(lags1[0])
    genotype_n_unsep1 = [genotype_n_sep1[i] + genotype_n_sep1[i + 1] for i in range(0, len(genotype_n_sep1), 2)]
    for s, species in enumerate(flask):
        for i, genotype in enumerate(species.genotypes):
            final_sub.append((seed, culture, rep, cycle, 1, species.name, genotype.name, genotype_n_sep1[2*i + (0, 2*nE)[s]], genotype_n_sep1[2*i + (1, 2*nE + 1)[s]], genotype_n_unsep1[i + (0, nE)[s]], genotype.lag, Ta, sol1[-1][0], sol1[-1][1], sol1[-1][2]))

    # mutation
    #genotype_n_sep_mut = generate_mutants(flask, copy.deepcopy(genotype_n_sep1), len(flask[0].genotypes), tuple([flask[s].mu for s, spec in enumerate(flask)]), mutation_function, cycle) #

    # init_cond2
    init_cond = list(init_R)
    for i in range(0, len(genotype_n_sep1), 2):
        init_cond.append(genotype_n_sep1[i])
        init_cond.append(genotype_n_sep1[i + 1])
    lags2 = [[i.lag for i in j.genotypes] for j in flask]

    # phase 2
    sol2 = run_phase(alpha, init_cond, lags2, t_grow, 2, rs=rs) ##

    # collect 2
    genotype_n_sep2 = list(sol2[-1, 3:])
    #genotype_n_unsep2 = [genotype_n_sep2[i] + genotype_n_sep2[i+1] for i in range(0, len(genotype_n_sep2), 2)]

    # mutate post
    genotype_n_sep_mut_pre = generate_mutants(flask, copy.deepcopy(genotype_n_sep2), len(flask[0].genotypes),
                                          tuple([flask[s].mu for s, spec in enumerate(flask)]), mutation_function,
                                          cycle)

    genotype_n_sep_mut = [[0, i][int(i > 0)] for i in genotype_n_sep_mut_pre] #<stdin>:1: DeprecationWarning: In future, it will be an error for 'np.bool_' scalars to be interpreted as an index IF NOT USING INT()

    genotype_n_unsep_mut = [genotype_n_sep_mut[i] + genotype_n_sep_mut[i+1] for i in range(0, len(genotype_n_sep_mut), 2)]

    # register final counts of genotypes in flask
    nE = len(flask[0].genotypes)
    # E
    for i, n in enumerate(genotype_n_unsep_mut[:nE]):   ##
        flask[0].genotypes[i].n = n / 2  # divide by 2 for transfer
    # S
    for j, n in enumerate(genotype_n_unsep_mut[nE:]): ## should not run on mono
        flask[1].genotypes[j].n = n / 2  # divide by 2 for transfer

    inher_R = (sol2[-1][0] / 2, sol2[-1][1] / 2, sol2[-1][2] / 2)

    # append 2
    for s, species in enumerate(flask):
        for i, genotype in enumerate(species.genotypes):
            final_sub.append((seed, culture, rep, cycle, 2, species.name, genotype.name, genotype_n_sep_mut[2*i + (0, 2*nE)[s]], genotype_n_sep_mut[2*i + (1, 2*nE + 1)[s]], genotype_n_unsep_mut[i + (0, nE)[s]], genotype.lag, Ta, inher_R[0], inher_R[1], inher_R[2])) ##

    return final_sub, inher_R

# simulation
def run(seed, culture, reps, mu, cycles, init_R, init_n, init_lag, Ta, alpha, t_grow, mutation_func_type, max_lag_change, rs):
    globals()['seed'] = seed ##
    globals()['rs'] = rs##

    #file = open(f'hcb_sim_{culture}_{seed}_met{init_R[0]}_lac{init_R[1]}.csv', 'w') # write custom text to front
    file = open(f'hcb_sim_{culture}_{seed}_met{init_R[0]}_rs{rs}.csv', 'w')  # write custom text to front
    file.write(f'##culture:{culture}#seed:{seed}#rep:{reps}#mu:{mu}#cycles:{cycles}#init_R:{init_R}#init_n:{init_n}#init_lag:{init_lag}#Ta:{Ta}#alpha:{alpha}#mut_func:{mutation_func_type}#max_lag_change:{max_lag_change}\n')

    # make mutation function
    if mutation_func_type == "null":
        mutation_func = make_null_function(max_lag_change)

    #print(f"Culture type: {culture}")#
    final = [('seed', 'culture', 'rep', 'cycle', 'phase_end', 'species', 'genotype', 'nlag', 'ngrow', 'ntot', 'lag', 'Ta', 'M', 'L', 'A')]
    for rep in range(reps):
        #print(f"Rep {rep}")#
        # set up first species
        flask = Flask()
        flask.add_species(Species('Escherichia coli', mu[0]))
        flask[0].add_genotype(Genotype('E0c0', init_n[0], init_lag[0], ancestors = '0'))
        final.append((seed, culture, rep, 0, 0, 'Escherichia coli', 'E0c0', init_n[0], 0, init_n[0], init_lag[0], Ta, init_R[0], init_R[1], init_R[2]))

        if culture == "co":
            flask.add_species(Species('Salmonella enterica', mu[1]))
            flask[1].add_genotype(Genotype('S0c0', init_n[1], init_lag[1], ancestors = '0'))
            final.append((seed, culture, rep, 0, 0, 'Salmonella enterica', 'S0c0', init_n[1], 0, init_n[1], init_lag[1], Ta, init_R[0], init_R[1], init_R[2]))

        inher_R = (0, 0, 0)
        # run simulation
        for cycle in range(cycles):
            #print(f"Cycle {cycle}")#
            final_sub, inher_R = run_one_simulation(seed, culture, flask, init_R, inher_R, Ta, alpha, t_grow, rep, cycle, mutation_func)
            for row in final_sub:
                final.append(row)

    final_pd = pd.DataFrame(final[1:], columns=list(final[0]))
    #with open(f'hcb_sim_{culture}_{seed}.csv', 'a') as file: # write custom text to front
    final_pd.to_csv(file, header=True, index=False, mode="a")
    #final_pd.to_csv(f'{culture}_seed{seed}_rep{reps}_mu{mu}_cycles{cycles}_init_R{init_R}_init_n{init_n}_init_lag{init_lag}_Ta{Ta}_alpha{alpha}_{mutation_func_type}{max_lag_change}.csv', index=False)

    #print("Finished")

#run(166, "co", 5, (0.0003, 0.0003), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), 0.5)
#run(166, "mono", 5, (0.0003, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))

#begin = time.perf_counter()  #
#run(499, "mono", 10, (0.0005, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), 0.5)
#print(f"{time.perf_counter() - begin}")  #
# resource, cycles, t_grow, first cycle