import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colormap

from scipy.integrate import odeint
#from tolerance_odes import odes
import seaborn as sns

import random
import copy
import math
import pdb#

## mutual extinction

def odes(x, t, alpha, n, tau_lag, frid=False): #
    '''
    x:[M, A, L, Egi, Eli, Sgi, Sli] initial conditions where Eg and El alternate indexes for i strains, then are followed by Sg and Sl alternating for j strains
    t:array of times
    alpha:[alphaE, alphaS] where alphaE: additional death rate due to antibiotic, as a proportion of max growth rate r for E. coli, alphaS: "" for S.enterica
    tau_lag:[tau_lagE, tau_lagS] where tau_lagE:list of tau_lags for each strain of E. coli, tau_lagS: "" for S. enterica

    n:[nE, nS] where nE:int strains of E. coli, nS: "" for S. enterica
    frid: when True, dMdt & dAdt & dLdt = 0 (for use with Fridman-style analysis)
    '''
    # read in parameters
    alphaE = alpha[0]
    alphaS = alpha[1]

    tau_lagE = tau_lag[0]
    tau_lagS = tau_lag[1]

    nE = n[0]
    nS = n[1]

    # R
    M = x[0]
    A = x[1]
    L = x[2]

    # half-saturation constants
    K_M = 1  #
    K_A = 1  #
    K_L = 1  #

    # resource decay constants
    kM = 5e-9  #
    kA = 5e-9  #
    kL = 5e-9  #

    # E
    for i in range(1, nE + 1):
        # growth constants
        locals()[f'rE{i}'] = 1  #
        locals()[f'kE{i}'] = 5e-9  #
        locals()[f'tau_lagE{i}'] = tau_lagE[i - 1]

        # resource constants
        locals()[f'cM{i}'] = 0.1  #
        locals()[f'pA{i}'] = 1.01  #
        locals()[f'cL{i}'] = 1.0  #

        # solutions -- starting with Eg1 on x[3] and El1 on x[4], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[1 + 2 * i]
        locals()[f'El{i}'] = x[2 + 2 * i]

        # differential equations
        locals()[f'dEg{i}dt'] = (1 - alphaE) * locals()[f'rE{i}'] * locals()[f'Eg{i}'] * (M / (M + K_M)) * (
                    L / (L + K_L)) - locals()[f'kE{i}'] * locals()[f'Eg{i}'] + locals()[f'El{i}'] / locals()[
                                    f'tau_lagE{i}']
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}'] / locals()[f'tau_lagE{i}']

    # S
    for j in range(1, nS + 1):
        # constants
        locals()[f'rS{j}'] = 0.5  #
        locals()[f'kS{j}'] = 5e-9  #
        locals()[f'tau_lagS{j}'] = tau_lagS[j - 1]

        # resource constants
        locals()[f'cA{j}'] = 1.0  #
        locals()[f'pM{j}'] = 1.56  #

        # solutions -- starting with Sg1 after the last El, Sg and Sl occupy alternating indexes for each strain
        locals()[f'Sg{j}'] = x[(1 + 2 * nE) + 2 * j]
        locals()[f'Sl{j}'] = x[(2 + 2 * nE) + 2 * j]

        # differential equations
        locals()[f'dSg{j}dt'] = (1 - alphaS) * locals()[f'rS{j}'] * locals()[f'Sg{j}'] * (M / (M + K_M)) * (
                    L / (L + K_L)) - locals()[f'kS{j}'] * locals()[f'Sg{j}'] + locals()[f'Sl{j}'] / locals()[
                                    f'tau_lagS{j}']
        locals()[f'dSl{j}dt'] = -locals()[f'Sl{j}'] / locals()[f'tau_lagS{j}']

    # M
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}'] * locals()[f'Eg{i}'] * (M / (M + K_M)) * (L / (L + K_L))
    sigma_pM = 0
    for j in range(1, nS + 1):
        sigma_pM += locals()[f'pM{j}'] * locals()[f'rS{j}'] * locals()[f'Sg{j}'] * (A / (A + K_A))

    dMdt = (-sigma_cM + sigma_pM - kM * M) * ([1, 0][frid])

    # A
    sigma_pA = 0
    for i in range(1, nE + 1):
        sigma_pA += locals()[f'pA{i}'] * locals()[f'rE{i}'] * locals()[f'Eg{i}'] * (M / (M + K_M)) * (L / (L + K_L))
    sigma_cA = 0
    for j in range(1, nS + 1):
        sigma_cA += locals()[f'cA{j}'] * locals()[f'Sg{j}'] * (A / (A + K_A))

    dAdt = (sigma_pA - sigma_cA - kA * A) * ([1, 0][frid])

    # L
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}'] * locals()[f'Eg{i}'] * (M / (M + K_M)) * (L / (L + K_L))

    dLdt = (-sigma_cL - kL * L) * ([1, 0][frid])

    to_return = [dMdt, dAdt, dLdt]
    for i in range(1, nE + 1):
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])
    for j in range(1, nS + 1):
        to_return.append(locals()[f'dSg{j}dt'])
        to_return.append(locals()[f'dSl{j}dt'])

    return to_return

### unedited, very messy

seed = random.randrange(1000)
random.seed(seed)
print('seed:', seed)

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

def run_tolerance(init_cond, lags, Ta, names_info): # gens for context
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors) #grabbed from stackoverflow

    # n
    nE = sum([i[0] == 'E' for i in names_info])
    nS = len(names_info) - nE
    n = [nE, nS]

    # lags undone
    lags_undo = [k for s in lags for k in s] # i still dont rly get why this works

    t_interval1 = np.linspace(0, Ta, 1000)
    sol1 = odeint(odes, init_cond, t_interval1, args=([2, 2], n, lags, False)) # false for now
    if globals()['genx'] == 0:
        j = 0
        for i in range(0, int(sol1[:, 3:].shape[1]), 2):
            plt.plot(t_interval1, np.add(sol1[:, i + 3], sol1[:, i + 4]), label= f'{round_half_up(lags_undo[j], 3)}, {names_info[j]}') # total pop
            j += 1
        plt.legend(loc='center right')
        plt.show()

    t_interval2 = np.linspace(0, 20, 1000) # arbitrary
    sol2 = odeint(odes, sol1[-1, :], t_interval2, args=([0, 0], n, lags, False))

    if globals()['genx'] == 0:
        j = 0
        for i in range(0, int(sol2[:, 3:].shape[1]), 2):
            plt.plot(t_interval2, np.add(sol2[:, i + 3], sol2[:, i + 4]), label= f'{round_half_up(lags_undo[j], 3)}, {names_info[j]}') # total pop
            j +=1
        plt.legend(loc='center right')
        plt.show()
    return sol2

def run_one_simulation(flask, init_R, Ta, rep, gen, mutation_function):
    final_sub = []
    for genotype in flask[0].genotypes:
        final_sub.append((rep, gen, flask[0].name, genotype.name, genotype.n, genotype.lag, Ta))
    for genotype in flask[1].genotypes:
        final_sub.append((rep, gen, flask[1].name, genotype.name, genotype.n, genotype.lag, Ta))

    # run phases 1, 2
    init_cond = init_R
    for genotype in flask[0].genotypes:
        init_cond.append(genotype.n)
        init_cond.append(0) # for each genotype, El = n and Eg = 0

    # S
    for genotype in flask[1].genotypes:
        init_cond.append(genotype.n)
        init_cond.append(0)

    lags = [[i.lag for i in flask[0].genotypes], [i.lag for i in flask[1].genotypes]]
    names_info = [i.name for i in flask[0].genotypes] + [i.name for i in flask[1].genotypes]
    sol = run_tolerance(init_cond, lags, Ta, names_info)

    # store per-gen data for graphing
    #per_gen_data_sub = (gen, tuple([i.name for i in flask[0].genotypes]), tuple(lags), sol)

    # counts after phases
    genotype_n_sep = [sol[-1, 3:len(flask[0].genotypes)*2 + 3], sol[-1, len(flask[0].genotypes)*2 + 3:]]#
    #pdb.set_trace()

    for flsk in range(2):#

        genotype_n_unsep = [ genotype_n_sep[flsk][i] + genotype_n_sep[flsk][i+1] for i in range(0, len(genotype_n_sep[flsk]), 2) ]
        genotype_freq = [ i/sum(genotype_n_unsep) for i in genotype_n_unsep ] ## faster to use numpy or sth?

        # mutation
        # Adamowicz-based
        chance = [random.uniform(0, 1) for i in range(round_half_up(sum(genotype_n_unsep)))]
        chance_tf = [i < flask[flsk].u for i in chance]
        mutant_n = sum(chance_tf)
        print(mutant_n)#

        if mutant_n != 0:
            #print('mutant_n', mutant_n)  # make sure mutant frequency is reasonable - 1's and 0's
            mutants = random.choices(range(len(genotype_n_unsep)), weights=genotype_freq, k=mutant_n)

            st_ct = len(flask[flsk].genotypes)
            for i, anc_i in enumerate(mutants):
                ancestor = flask[flsk].genotypes[anc_i]

                genotype_n_unsep[anc_i] -= 1

                pref = ['E', 'S'][flsk]
                flask[flsk].add_genotype(Genotype(f'{pref}{st_ct + i}g{gen}', n = 1, lag = max([0, ancestor.lag + mutation_function()]), ancestors = ancestor.name + ' ' + ancestor.ancestors))

        # register final counts of genotypes in Flask flask
        for i in range(len(genotype_n_unsep)):
            flask[flsk].genotypes[i].n = genotype_n_unsep[i]

    return final_sub

# keep for now -----
def tuple_list_to_pd_dataframe(tuple_list):
    dic = {}
    for ind in range(len(tuple_list[0])):
        dic[tuple_list[0][ind]] = [i[ind] for i in tuple_list[1:]]
        
    return pd.DataFrame(dic)
# keep for now -----

### simulation parameters
reps = 1
u = 0.001 # mutation rate
gens = 50

init_R = [0.01, 0.01, 100] # starting [M, A, L] of each new growth flask ## should it be just lactose?
Ta = 3 # length of anibiotic treatment
max_lag_change = 1.1 # max mutation-induced lag change ## orig. antibiotic_change_per_well * 1.1

### make mutation function
# null_function grabs MIC mutations from a uniform distribution
null_function = make_null_function(max_lag_change)

### run simulation
final = [('rep', 'gen', 'species', 'genotype', 'n', 'lag', 'Ta')]
per_gen_data = []
for rep in range(reps):

    # set up first species
    flask = Flask()
    flask.add_species(Species('Escherichia coli', u))
    ## lag actual value does not affect - 3
    flask[0].add_genotype(Genotype('E0g0', n = 1, lag = 1, ancestors = '0')) ## lag cannot be 0 to avoid divide by zero

    # set up second species
    flask.add_species(Species('Salmonella enterica', u)) # same u as E. coli for now
    flask[1].add_genotype(Genotype('S0g0', n=1, lag=1, ancestors='0'))

    # run simulation
    for gen in range(gens):
        globals()['genx'] = gen
        final_sub = run_one_simulation(flask, copy.deepcopy(init_R), Ta, rep, gen, null_function)
        for row in final_sub:
            final.append(row)
        #per_gen_data.append(per_gen_data_sub)

final_pd = tuple_list_to_pd_dataframe(final)
final_pd.to_csv('final_co.csv', index=False)

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