import pdb
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tolerance_odes_heatmap_rev import odes

import random
import copy
import math

# make seed visible for re-generating purposes
# globals()['seed'] = random.randrange(1000)
# random.seed(seed)
# print('seed:', seed)

# rounding function: rounds .5 up to next int
# def round_half_up(n, decimals=0):
#     multiplier = 10 ** decimals
#     rounded = math.floor(n * multiplier + 0.5) / multiplier
#     return int(rounded)

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

def run_phase(alpha, init_cond, lags, t, phase, inc=1000, frid=False): #
    alpha_this = tuple([[a,0][phase-1] for a in alpha])
    t_interval = np.linspace(0, t, inc)
    sol = odeint(odes, init_cond, t_interval, args=(alpha_this, lags, frid)) ##

    return sol

def generate_mutants(flask, n, mu, mutation_function, cycle):
    nE = len(flask[0].genotypes)
    if len(flask) > 1:
        nS = len(flask[1].genotypes)
    else:
        nS = 0

    n_mut_grow = [[], []]
    n_mut_lag = [[], []]
    for s, species in enumerate(flask):
        # pdb.set_trace()
        n_growing = n[s*(2*nE):(nE, 2*nE + nS)[s]] # pick sol growing

        u = mu[s]
        sum = np.sum(n_growing)
        n_freq = n_growing / sum

        # Adamowicz-based
        chance_tf = np.random.rand(int(sum))
        mutant_n = len(chance_tf[chance_tf < u])

        if mutant_n > 0:
            mutants = random.choices(range(len(n_growing)), weights=n_freq, k=mutant_n)
            genotype_ct = len(flask[s].genotypes)

            for j, anc_i in enumerate(mutants):
                ancestor = flask[s].genotypes[anc_i]

                n[anc_i + (0, 2*nE)[s]] -= 1 # tuple is for correct indices in n, which is not separated by species
                n_mut_grow[s].append(1)
                n_mut_lag[s].append(0)
                # n = np.insert(n, (nE, 2*nE + nS)[s], 1)
                # n = np.insert(n, 2*nE + s*(2*nS) + j + 1, 0) # insert after latest lagging
                flask[s].add_genotype(Genotype(f"{('E', 'S')[s]}{genotype_ct + j}c{cycle}", n=1, lag=max([0, ancestor.lag + mutation_function(s)]),
                                               ancestors=ancestor.name + ' ' + ancestor.ancestors))

    return n_mut_grow, n_mut_lag

def run_one_simulation(seed, culture, flask, init_R, inher_R, Ta, alpha, t_grow, rep, cycle, mutation_function):
    final_sub = []

    # init_cond1
    init_cond = [init_R[0] + inher_R[0], init_R[1] + inher_R[1], init_R[2] + inher_R[2]]

    for species in flask:
        # should not reset that which is already appended
        Gs = []
        Ls = []
        for genotype in species.genotypes:
            Gs.append(0)
            Ls.append(genotype.n)
        init_cond += Gs + Ls
    init_cond = np.array(init_cond)#
    lags1 = [[i.lag for i in j.genotypes] for j in flask]

    # phase 1
    sol = run_phase(alpha, init_cond, lags1, Ta, 1) ##

    # collect 1
    sol1 = sol[-1]

    # append 1
    nE = len(lags1[0])
    for s, species in enumerate(flask):
        N = len(lags1[s])
        N_growing = sol1[3 + s*2*nE:3 + nE + s*(nE + N)]
        N_lagging = sol1[3 + nE + s*(nE + N):]
        for i, genotype in enumerate(species.genotypes):
            final_sub.append((culture, rep, cycle, 1, species.name, genotype.name, N_growing[i], N_lagging[i], genotype.lag, sol1[0], sol1[1], sol1[2]))

    # init_cond2
    init_cond = list(init_R) + list(sol1[3:]) ## gross
    lags2 = [[i.lag for i in j.genotypes] for j in flask]

    # phase 2
    sol = run_phase(alpha, init_cond, lags2, t_grow, 2) ##

    # collect 2
    sol2 = sol[-1]

    # mutate post
    n = sol2[3:]
    nE = len(flask[0].genotypes)
    n_mut_grow, n_mut_lag = generate_mutants(flask, n,
                                          tuple([flask[s].mu for s, spec in enumerate(flask)]), mutation_function,
                                          cycle)

    # append 2 and update flask at same time
    for s, species in enumerate(flask):
        N = [n[0:2*nE], n[2*nE:]][s]
        half = (int)(len(N) / 2)
        # pdb.set_trace()

        N_growing = list(N[:half])
        N_growing.extend(n_mut_grow[s])
        N_lagging = list(N[half:])
        N_lagging.extend(n_mut_lag[s])
        
        N_tot = np.array(N_growing) + np.array(N_lagging)

        for i, ct in enumerate(N_tot):
            #pdb.set_trace()
            genotype = species.genotypes[i]     ## Ntot more than flask
            genotype.n = ct / 2

            final_sub.append((culture, rep, cycle, 2, species.name, genotype.name, N_growing[i], N_lagging[i], genotype.lag, sol2[0], sol2[1], sol2[2]))

    inher_R = (sol2[0] / 2, sol2[1] / 2, sol2[2] / 2)

    return final_sub, inher_R

# simulation
def run(seed, culture, reps, mu, cycles, init_R, init_n, init_lag, Ta, alpha, t_grow, mutation_func_type, max_lag_change):
    globals()['seed'] = seed ##

    #file = open(f'hcb_sim_{culture}_{seed}_met{init_R[0]}_lac{init_R[1]}.csv', 'w') # write custom text to front
    file = open(f'hcb_sim_{culture}_{seed}_met{init_R[0]}.csv', 'w')  # write custom text to front
    file.write(f'##culture:{culture}#seed:{seed}#rep:{reps}#mu:{mu}#cycles:{cycles}#init_R:{init_R}#init_n:{init_n}#init_lag:{init_lag}#Ta:{Ta}#alpha:{alpha}#mut_func:{mutation_func_type}#max_lag_change:{max_lag_change}\n')

    # make mutation function
    if mutation_func_type == "null":
        mutation_func = make_null_function(max_lag_change)

    #print(f"Culture type: {culture}")#
    final = [('culture', 'rep', 'cycle', 'phase_end', 'species', 'genotype', 'ngrow', 'nlag', 'lag', 'M', 'L', 'A')]
    for rep in range(reps):
        #print(f"Rep {rep}")#
        # set up first species
        flask = Flask()
        flask.add_species(Species('Escherichia coli', mu[0]))
        flask[0].add_genotype(Genotype('E0c0', init_n[0], init_lag[0], ancestors = '0'))
        final.append((culture, rep, 0, 0, 'Escherichia coli', 'E0c0', 0, init_n[0], init_lag[0], init_R[0], init_R[1], init_R[2]))

        if culture == "co":
            flask.add_species(Species('Salmonella enterica', mu[1]))
            flask[1].add_genotype(Genotype('S0c0', init_n[1], init_lag[1], ancestors = '0'))
            final.append((culture, rep, 0, 0, 'Salmonella enterica', 'S0c0', 0, init_n[1], init_lag[1], init_R[0], init_R[1], init_R[2]))

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


#run(166, "mono", 5, (0.0003, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
run(random.randrange(1000), "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

# begin_tot = time.perf_counter()  #
#run(499, "mono", 10, (0.0005, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), 0.5)
# print(f"{time.perf_counter() - begin_tot}")  #
#run(499, "co", 5, (0.0003, 0.0003), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), 0.5)