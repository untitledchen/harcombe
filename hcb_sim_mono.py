##### change only calculate fitness to run odes and return final proportions as w
##### remove olds like MIC and s? s = 1 or r or whatever
##### question on the original code == overwrite? 289 to 103 model
##### question what does the t in odes() do? ODES folder

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#
from scipy.integrate import odeint

import random
import copy
import math
import pdb#

seed = random.randrange(1000)
random.seed(seed)
print('seed:', seed)

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

##### mutation_functions.r
def make_null_function(chance_of_MIC_mutation, max_MIC_change, sd_s_change):
    def null_function():
        if random.uniform(0,1) < chance_of_MIC_mutation:
            MIC = random.uniform(0,1) * max_MIC_change
            s = 0
        else:
            MIC = 0
            s = random.normalvariate(0, sd_s_change)
        return (s, MIC)
    return null_function
            
##### model_functions2.r
class Genotype():
    def __init__(self, name, n, s=1, MIC=0, ancestors='0'):
        self.name = name
        self.n = n
        self.s = s
        self.MIC = MIC
        self.ancestors = ancestors

    def __str__(self):
        return "{'name': " + str(self.name) + ", 'n': " + str(self.n) + ", 's': " + str(self.s) + ", 'MIC': " + str(self.MIC) + ", 'ancestors': " + str(self.ancestors) +"}"

class Species():
    def __init__(self, name, N, u):
        self.name = name
        self.N = N
        self.u = u
        self.genotypes = []

    def __str__(self):
        return "{'genotypes': " + str([i.name for i in self.genotypes]) + ", 'N': " + str(self.N) + ", 'u': " + str(self.u) + "}"

    def add_genotype(self, n, s=1, MIC=0, ancestors='0', genotype_name=None, genotype_num=1):
        if genotype_num == 1:
            if genotype_name == None:
                genotype_name = len(self.genotypes)
                
            if genotype_name in [i.name for i in self.genotypes]:
                ind = [i.name for i in self.genotypes].index(genotype_name)
                #self.genotypes[ind].n += n
                self.genotypes[ind] = Genotype(genotype_name, n, s, MIC, ancestors)
            else:
                self.genotypes.append(Genotype(genotype_name, n, s, MIC, ancestors))
        else:
            if genotype_name == None:
                genotype_name = [(len(self.genotypes) + i) for i in range(genotype_num)]

            for i in range(genotype_num):
                if genotype_num[i] in [i.name for i in self.genotypes]:
                    ind = [i.name for i in self.genotypes].index(genotype_num[i])
                    #self.genotypes[ind].n += n[i]
                    self.genotypes[ind] = Genotype(genotype_name[i], n[i], s[i], MIC[i], ancestors[i])
                else:
                    self.genotypes.append(Genotype(genotype_name[i], n[i], s[i], MIC[i], ancestors[i]))

class Well(list):
    def __init__(self, name):
        list.__init__(self)
        self.name = name

    def __str__(self):
        return str([i.name for i in self])

    def add_species(self, N, u, species_name=None, species_num=1):
        if species_num == 1:
            if species_name == None:
                species_name = len(self)
            self.append(Species(species_name, N, u))
        else:
            if species_name == None:
                species_name = [(len(self) + i) for i in range(species_num)]
            for i in range(species_num):
                self.append(Species(species_name[i], N[i], u[i]))
    
class Season(list):
    def __init__(self):
        list.__init__(self)

    def __str__(self):
        return str([i.name for i in self])

    def add_well(self, well_name=None, well_num=1):
        if well_num == 1:
            if well_name == None:
                well_name = len(self)
            self.append(Well(well_name))
        else:
            if well_name == None:
                well_name = [(len(self) + i) for i in range(well_num)]
            for i in range(well_num):
                self.append(Well(well_name[i]))

def make_first_season(wells, n_species, N, u, genotype_prefixes):
    first_season = Season()

    first_season.add_well(well_num=wells)
    for well in range(wells):
        for k in range(n_species):
            first_season[well].add_species(N, u)
            first_season[well][k].add_genotype(N, genotype_name=genotype_prefixes[well]+'_0')

    return first_season

def get_n(Species):
    return [i.n for i in Species.genotypes]
'''
def get_s(Species):
    return [i.s for i in Species.genotypes]

def get_MIC(Species):
    return [i.MIC for i in Species.genotypes]
'''

## road work ahead
def calculate_fitness(s, n, MIC, antibiotic):
    MIC_tf = [i < antibiotic for i in MIC]
    new_s = [[i,0][true] for i,true in zip(s, MIC_tf)]
    
    if all([i == 0 for i in new_s]):
        print('All individuals have zero fitness')
        return -1

    w_pre = [i/sum(new_s)*j for i, j in zip(new_s, n)]
    w_next = [i/sum(w_pre) for i in w_pre]
    w = pd.Series(w_next).cumsum()
    return list(w)

def calculate_fitness(init_conc, species, Ta):
    '''
    init_conc:list [init_M, init_L]
    species:Species
    Ta:int antibiotic duration in hr

    return relative proportions of genotypes after phase 1 and 2, as cumulative sum
    '''

def run_phase(odes, init_cond, t_interval, init_lag, phase, frid=False):
    if phase == 1:
        alpha = 2
    elif phase == 2:
        alpha = 0
    
    sol = odeint(odes, init_cond, t_interval, args=(alpha, init_lag, frid)) ## when alpha = -2, effective growth rate = -r
    para = [t_interval, init_lag] # removed strain
    return sol, para

def regular_mono(init_res, init_freq, init_lag, Ta, frid=False):
    '''
    init_conc:list [init_M, init_L]
    init_freq:list of frequencies of each genotype
    init_lag:list of lag times of each genotype
    Ta: length of killing phase
    frid: when True, dMdt and dLdt = 0 (for determining effective starting population)
    '''
    # init_cond
    pops = [0]*len(init_freq) * 2
    pops[1::2] = init_freq
    
    init_cond = [init_res] + pops
    
    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, Ta*2, 1000) ## arbitrary

    # run phases
    sol1, para1 = run_phase(odes, init_cond, t_interval1, init_lag, phase=1)
    sol2, para2 = run_phase(odes, sol1[-1, :], t_interval2, init_lag, phase=2, frid=frid)

    # convert sols to w
    totals = [sol2[i] + sol2[i + 1] for i in range(2, len(sol2), 2)] # confirmed
    totals_proportion = [i/sum(totals) for i in totals]
    w = pd.Series(totals_proportion).cumsum() ##

    return list(w)

## road work ahead

def birth_offspring(N, w):
    chance = [random.uniform(0,1) for i in range(N)]
    return [[i < j for j in w].index(True) for i in chance]

def get_mutants(Species):
    chance = [random.uniform(0,1) for i in range(Species.N)]
    chance_tf = [i < Species.u for i in chance]
    return [ind for ind, true in enumerate(chance_tf) if true]

def make_mutants(Species, ancestor_genotype_names, new_genotype_names, mutant_function):
    for i in range(len(ancestor_genotype_names)):
        ancestor_name = ancestor_genotype_names[i][1]
        ancestor = Species.genotypes[ancestor_genotype_names[i][0]]
        
        Species.add_genotype(genotype_name=new_genotype_names[i], n=1, s=max([0, ancestor.s+mutant_function()[0]]), MIC=max([0, ancestor.MIC+mutant_function()[1]]), ancestors= ancestor_name+' '+ancestor.ancestors)
    return Species

def set_n(Species, genotype_ns):
    if sum(genotype_ns) != Species.N: #
        print('Sum of n != N defined in Species')
        return -1
    
    for i in range(len(Species.genotypes)):
        Species.genotypes[i].n = genotype_ns[i]
    return Species

def run_one_simulation(species, gens, antibiotic, mutant_func, interdependent, current_prefix):
    n_species = len(species)

    antibiotic = [antibiotic]*gens

    results = [[('gen', 'genotypes', 'freq', 's', 'MIC', 'ancestors')] for i in range(n_species)]
    
    if interdependent:
        if any([len(i.genotypes) == 0 for i in species]):
            for k in range(n_species):
                results[k].append((-1, '-1', -1, -1, -1, '-1'))
            return results

    for k in range(n_species):
        for i in species[k].genotypes:
            results[k].append((0, i.name, i.n / species[k].N, i.s, i.MIC, i.ancestors))

    for gen in range(1, gens+1):
        globals()['genx'] = gen#
        
        w = [calculate_fitness(get_s(k), get_n(k), get_MIC(k), antibiotic[gen-1]) for k in species]
        
        if interdependent:
            if any([i == -1 for i in w]):
                for k in range(n_species):
                    results[k].append((-1, '-1', -1, -1, -1, '-1'))
                return results

        pop_sizes = [i.N for i in species]
        offspring = [birth_offspring(pop_sizes[k], w[k]) for k in range(n_species)]

        mutants = [get_mutants(k) for k in species]

        for k in range(n_species):
            if len(mutants[k]) > 0:
                genotype_names = [i.name for i in species[k].genotypes]
                max_genotype_name = max([int(i.split('_')[-1]) for i in genotype_names])

                new_genotype_names = [f'{current_prefix}_{max_genotype_name + i}' for i in range(1, len(mutants[k])+1)]
                ancestor_genotype_names = [(offspring[k][i], species[k].genotypes[offspring[k][i]].name) for i in mutants[k]]

                species[k] = make_mutants(species[k], ancestor_genotype_names, new_genotype_names, mutant_function)
        
                for i in range(len(new_genotype_names)):
                    offspring[k][mutants[k][i]] = len(genotype_names) + i

        genotype_ns = [[offspring[k].count(i) for i in range(len(species[k].genotypes))] for k in range(n_species)]
        species = [set_n(species[k], genotype_ns[k]) for k in range(n_species)]
       
        for k in range(n_species):
            for i in species[k].genotypes:
                results[k].append((gen, i.name, i.n / species[k].N, i.s, i.MIC, i.ancestors))

    return results

def start_next_season(results, wells, n_species, N, u):
    next_season = Season()
    next_season.add_well(well_num=wells)
    
    for well in range(wells):
        if well == 0:
            for k in range(n_species):
                max_gen = results[well][k][-1][0]
                curr_result = [i for i in results[well][k] if (i[0] == max_gen and i[2] > 0)]

                genotype_names = [(ind, val[1]) for ind, val in enumerate(curr_result)]
                genotype_ns = [i[2]*N for i in curr_result]
                
                genotypes_pre = [[genotype_names[i][1]]*round_half_up(genotype_ns[i]) for i in range(len(genotype_names))]
                genotypes = [element for sublist in genotypes_pre for element in sublist]

                new_genotype_names = genotype_names
                new_genotype_ns = genotype_ns

                next_season[well].add_species(N, u)
                for i in range(len(new_genotype_names)):
                    ind = new_genotype_names[i][0]
                    next_season[well][k].add_genotype(new_genotype_ns[i], curr_result[ind][3], curr_result[ind][4], curr_result[ind][5], genotype_name=new_genotype_names[i][1])
        else:
            for k in range(n_species):
                next_season[well].add_species(N, u)
                draw_wells = [results[well-1][k], results[well][k]]
                extinction = [any([(i[2] == -1) for i in draw_wells[0]]), any([(i[2] == -1) for i in draw_wells[1]])]
                if all(extinction):
                    cells_per_well = [0, 0]
                elif extinction[0] and not extinction[1]:
                    cells_per_well = [0, N]
                elif extinction[1] and not extinction[0]:
                    cells_per_well = [N, 0]
                else:
                    cells_per_well = [N/2, N/2]
                for draw_well in range(2):
                    if cells_per_well[draw_well] == 0:
                        continue
                    max_gen = draw_wells[draw_well][-1][0]
                    curr_result = [i for i in draw_wells[draw_well] if (i[0] == max_gen and i[2] > 0)]
 
                    genotype_names = [(ind, val[1]) for ind, val in enumerate(curr_result)]
                    genotype_ns = [i[2]*N for i in curr_result]
                    
                    genotypes_pre = [[genotype_names[i]]*round_half_up(genotype_ns[i]) for i in range(len(genotype_names))]
                    genotypes = [element for sublist in genotypes_pre for element in sublist]
                    genotypes = random.sample(genotypes, round_half_up(cells_per_well[draw_well]))

                    new_genotype_names = pd.Series(genotypes).unique()
                    new_genotype_names = list(new_genotype_names)
                    new_genotype_names.sort()
                    new_genotype_ns = [genotypes.count(i) for i in new_genotype_names]

                    for i in range(len(new_genotype_names)):
                        ind = new_genotype_names[i][0]
                        next_season[well][k].add_genotype(new_genotype_ns[i], curr_result[ind][3], curr_result[ind][4], curr_result[ind][5], genotype_name=new_genotype_names[i][1])

    return next_season

def summarize_results(results, wells, n_species, season, u, rep, gens, mutant_function):
    final = []
    for well in range(wells):
        for k in range(n_species):
            curr_result = results[well][k]
            if any([i[0] == -1 for i in curr_result]):
                final.append((well, k, season, False, 0, np.nan, np.nan, np.nan, n_species, u, rep, gens, mutant_function))
                continue
            max_gen = curr_result[-1][0]
            curr_result = [i for i in curr_result if (i[0] == max_gen and i[2] > 0)]

            final.append((well, k, season, True, len(curr_result), sum([i[3]*i[2] for i in curr_result]), sum([i[4]*i[2] for i in curr_result]), sum([(len(i[5].split(' '))-1)*i[2] for i in curr_result]), n_species, u, rep, gens, mutant_function))## n_genotypes originally = len(curr_result['genotypes'])
    return final

def tuple_list_to_pd_dataframe(tuple_list):
    dic = {}
    for ind in range(len(tuple_list[0])):
        dic[tuple_list[0][ind]] = [i[ind] for i in tuple_list[1:]]
        
    return pd.DataFrame(dic)
        
##### adamowicz_et_al_evolution_model_code_example.r
### simulation parameters
reps = 5 # reps per treatment condition
N = 1000 # individuals per species
u = 0.001 # mutation rate
max_interdependent_species = 3
seasons = 20 # number of transfers
gens = 20 # gens per well per season
wells = 15
antibiotic_change_per_well = 1 # antibiotic concentration increase per well

### mutation function parameters
chance_of_MIC_mutation = 0.5 # probability that a mutation is a mutation in MIC
max_MIC_change = antibiotic_change_per_well * 1.1 # max mutation-induced MIC change
sd_s_change = 0.1 # if a mutation affects growth rate, the change is a normally-distributed variable centered on zero with this sd

### make mutation function
# null_function grabs MIC mutations from a uniform distribution
null_function = make_null_function(chance_of_MIC_mutation = chance_of_MIC_mutation, max_MIC_change = max_MIC_change, sd_s_change = sd_s_change)

### choose which mutation function to use
mutant_function = null_function
mutant_function_name = 'null'

### model
n_species_in_consortium = [i for i in range(1, max_interdependent_species+1)]
antibiotic = [i*antibiotic_change_per_well for i in range(wells)]
current_prefixes = [f'w{well}s0' for well in range(wells)]

all_data = pd.DataFrame()
final = [('well', 'species', 'season', 'alive', 'n_genotypes', 's', 'MIC', 'mutations', 'n_species', 'u', 'rep', 'gens', 'mutant_function')]
for rep in range(reps):
    for n_species in n_species_in_consortium:
        print('rep', rep)
        print('\tgens', gens)
        print('\tu', u)
        print('\tmutant func', mutant_function_name)
        print('\tn_species', n_species)
        next_season = make_first_season(wells, n_species, N, u, genotype_prefixes=current_prefixes)

        for season in range(seasons):
            print('season', season)#
            globals()['seasonx'] = season#
            results = [run_one_simulation(next_season[well], gens, antibiotic[well], mutant_function, True, current_prefixes[well]) for well in range(wells)]
            next_season = start_next_season(copy.deepcopy(results), wells, n_species, N, u)
            current_prefixes = [f'w{well}s{season}' for well in range(wells)]
            final = final + summarize_results(copy.deepcopy(results), wells, n_species, season, u, rep, gens, mutant_function_name)

all_data = tuple_list_to_pd_dataframe(final)
all_data.to_csv('all_data.csv', index=False)

###
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
