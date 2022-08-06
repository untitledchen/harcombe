import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import random
import copy
import math
import pdb#

seed = random.randrange(1000)
random.seed(seed)
print('seed:', seed)
#"{'name': " + str(self.name) + "'n': " + str(self.n) + "'s': " + str(self.s) + "'MIC': " + str(self.MIC) + "'ancestors': " + str(self.ancestors) +"}"
#"{'genotypes': " + str(self.genotypes) + "'N': " + str(self.N) + "'u': " + str(self.u) + "}"

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
    def __init__(self, curr_prefix, N):
        self.name = f'{curr_prefix}_0'
        self.n = N
        self.s = 1
        self.MIC = 0
        self.ancestors = '0'

    def __str__(self):
        return "{'name': " + str(self.name) + ", 'n': " + str(self.n) + ", 's': " + str(self.s) + ", 'MIC': " + str(self.MIC) + ", 'ancestors': " + str(self.ancestors) +"}"


class Species():
    def __init__(self, k, N, u, genotype):
        self.name = k
        self.N = N
        self.u = u
        self.genotypes = [genotype]

    def __str__(self):
        return "{'genotypes': " + str([i.name for i in self.genotypes]) + ", 'N': " + str(self.N) + ", 'u': " + str(self.u) + "}"

class Well():
    def __init__(self, well, n_species, N, u, curr_prefix):
        self.name = well
        self.n_species = []
        for k in range(n_species):
            self.n_species.append(Species(k, N, u, Genotype(curr_prefix, N)))

    def __str__(self):
        return str([i.name for i in self.n_species])
    
class Season():
    def __init__(self, wells, prefixes, n_species, N, u):
        self = []

        for well in range(wells):
            curr_prefix = prefixes[well]
            self.append(Well(well, n_species, N, u, curr_prefix))

    def __iter__(self):
            yield well

    def __str__(self):
        return str([i.name for i in self])
##override iter
def Species(N=1000, u=0.001, first_genotype=None):
    me = {'genotypes':list(), 'N':N, 'u':u}##

    if first_genotype != None:
        if {'name', 'n', 's', 'MIC', 'ancestors'}.issubset(set(first_genotype[0])):
            me['genotypes'][first_genotype['name']] = first_genotype
        else:
            print('All genotypes must contain the following fields: name (character), n (integer), s (number), MIC (number), ancestors (character)')
            return 0
        
    return me

def make_first_season(wells, n_species, current_prefixes, u=0.0025):
    next_season = {}
    for well in range(1, wells+1):
        current_prefix = current_prefixes[well-1]
        first_genotype = {'name': f'{current_prefix}_0', 'n': N, 's': 1, 'MIC': 0, 'ancestors': '0'}
        
        species = {}
        for k in range(1, n_species+1):
            species[f'{k}'] = Species(N, u, first_genotype.copy(deep=True))

        next_season[f'{well}'] = species

    return next_season

def get_mutants(Species):
    if type(Species) == SpeciesType:
        offspring = pd.Series([random.uniform(0,1) for i in range(Species['N'])])
        mutants = pd.Series([i for i in offspring[offspring < Species['u']].index], dtype=object)
        return mutants
    else:
        print('This method requires a SpeciesType object')
        return Species

def add_genotype(Species, new_genotype):
    if type(Species) == SpeciesType:
        Species['genotypes'][new_genotype['name']] = new_genotype
        return Species
    else:
        print('This method requires a SpeciesType object')
        return Species

def set_n(Species, n):
    if type(Species) == SpeciesType:
        if sum(n) != Species['N']:
            print('Sum of n != N defined in Species')
            return -1
    
        for i in range(len(Species['genotypes'])):
            Species['genotypes'][i]['n'] = n[i]
        return Species
    else:
        print('This method requires a SpeciesType object')
        return Species

def make_mutants(Species, ancestor_genotype_names, new_genotype_names, mutant_function):
    if type(Species) == SpeciesType:
        if len(ancestor_genotype_names) == 0:
            print('Must supply ancestor_genotype_names to make mutants')
            return -1
        for k in range(len(ancestor_genotype_names)):
            ancestor_name = ancestor_genotype_names[k]
            ancestor = Species['genotypes'][ancestor_name]
            mutant = pd.Series({'name':new_genotype_names[k], 'n':1, 's':max([0, float(ancestor['s']+ mutant_function()[0])]), 'MIC':max([0, float(ancestor['MIC']+ mutant_function()[1])]), 'ancestors':ancestor_name+' '+ancestor['ancestors']})
            
            Species = add_genotype(Species, mutant.copy(deep=True))
        return Species
    else:
        print('This method requires a SpeciesType object')
        return Species

def birth_offspring(N, w):
    chance = pd.Series([random.uniform(0,1) for i in range(N)])
    return pd.Series([w[x < w].index[0] for x in chance])

def start_next_season(results, wells, n_species, N, u=0.0025):
    next_season = pd.Series(dtype=object)
    for well in range(1, wells+1):
        next_season[f'{well}'] = pd.Series(dtype=object)
        if well == 1:
            for k in range(1, n_species+1):
                curr_result = results['1'][f'{k}'].loc[results['1'][f'{k}']['gen'] == max(results['1'][f'{k}']['gen'])].loc[results['1'][f'{k}']['freq'] > 0]
                #curr_result['ancestors'] = as.character

                genotype_names = [i for i in curr_result['genotypes'].values]
                genotype_ns = (N * curr_result['freq']).reset_index(drop=True)
                genotypes = pd.Series([[genotype_names[x]]*int(genotype_ns[x]) for x in range(len(genotype_names))]).explode()

                new_genotype_names = pd.Series(genotypes.unique())
                new_genotype_ns = pd.Series([sum(genotypes == x) for x in new_genotype_names])

                new_population = Species(N, u)
                for i in range(len(new_genotype_names)):
                    new_genotype = pd.Series({'name': new_genotype_names[i], 'n': new_genotype_ns[i], 's': curr_result['s'].loc[curr_result['genotypes'] == new_genotype_names[i]], 'MIC': curr_result['MIC'].loc[curr_result['genotypes'] == new_genotype_names[i]], 'ancestors': curr_result['ancestors'].loc[curr_result['genotypes'] == new_genotype_names[i]]})
                    new_population = add_genotype(new_population, new_genotype)
                next_season[f'{well}'][f'{k}'] = new_population
        else:
            for k in range(1, n_species+1):
                new_population = Species(N, u)
                draw_wells = pd.Series([0], dtype=object)
                draw_wells[0] = results[f'{well-1}'][f'{k}']
                draw_wells[1] = results[f'{well}'][f'{k}']
                extinction = [any(draw_wells[0]['freq'] == -1), any(draw_wells[1]['freq'] == -1)]
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
                    curr_result = draw_wells[draw_well].loc[draw_wells[draw_well]['gen'] == max(draw_wells[draw_well]['gen'])].loc[draw_wells[draw_well]['freq'] > 0]
                    #curr_result['ancestors'] = as.character
                    
                    genotype_names = [i for i in curr_result['genotypes'].values]
                    genotype_ns = (N * curr_result['freq']).reset_index(drop=True)
                    genotypes = pd.Series([[genotype_names[x]]*round(genotype_ns[x]) for x in range(len(genotype_names))]).explode() ## better idea than int()?
                    genotypes = genotypes.sample(int(cells_per_well[draw_well]))

                    new_genotype_names = pd.Series(genotypes.unique())
                    new_genotype_ns = pd.Series([sum(genotypes == x) for x in new_genotype_names])

                    for i in range(len(new_genotype_names)):
                        new_genotype = pd.Series({'name': new_genotype_names[i], 'n': new_genotype_ns[i], 's': curr_result['s'].loc[curr_result['genotypes'] == new_genotype_names[i]], 'MIC': curr_result['MIC'].loc[curr_result['genotypes'] == new_genotype_names[i]], 'ancestors': curr_result['ancestors'].loc[curr_result['genotypes'] == new_genotype_names[i]]})
                        new_population = add_genotype(new_population, new_genotype)
                next_season[f'{well}'][f'{k}'] = new_population
            
    return next_season

def get_n(Species):
    if type(Species) == SpeciesType:
        return [i['n'] for i in Species['genotypes']]
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_s(Species):
    if type(Species) == SpeciesType:
        return [i['s'] for i in Species['genotypes']]
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_MIC(Species):
    if type(Species) == SpeciesType:
        return [i['MIC'] for i in Species['genotypes']]
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_ancestors(Species):
    if type(Species) == SpeciesType:
        return [i['ancestors'] for i in Species['genotypes']]
    else:
        print('This method requires a SpeciesType object')
        return Species

def calculate_fitness(s, n, MIC, antibiotic):
    MIC_tf = [i < antibiotic for i in MIC]
    s = [[i,0][true] for i,true in zip(s, MIC_tf)]
    
    if all(i == 0 for i in s):
        print('All individuals have zero fitness')
        return -1
    ##
    w = [i/sum(s) for i in s]
    w = pd.Series([i/sum(s) for i in s]).cumsum()
    return list(w)

def run_one_simulation(species, gens, antibiotic, mutant_func, interdependent=True, current_prefix=''):
    n_species = len(species)
    pop_sizes = [i['N'] for i in species]

    antibiotic = [antibiotic]*gens # removed the len check

    results = {} #moved it up here
    if interdependent:
        if any([len(i['genotypes']) == 0 for i in species]):
            for k in range(1, n_species+1):
                results[f'{k}'] ={'gen':-1, 'genotypes':'-1', 'freq':-1, 's':-1, 'MIC':-1, 'ancestors':'-1'}
            return results

    for k in range(1, n_species+1):
        results[f'{k}'] = {'gen':[0], 'genotypes': [[i for i in species[f'{k}']['genotypes'].keys()]], 'freq': [[i/species[f'{k}']['N'] for i in get_n(species[f'{k}'])]], 's': [get_s(species[f'{k}'])], 'MIC': [get_MIC(species[f'{k}'])], 'ancestors': [get_ancestors(species[f'{k}']])} ##

    for gen in range(1, gens+1):
        w = [calculate_fitness(get_s(k), get_n(k), get_MIC(k), antibiotic[gen-1]) for k in species] ##
        if interdependent:
            if any(w == -1):
                for k in range(1, n_species+1):
                    results[f'{k}'] = {'gen':-1, 'genotypes':'-1', 'freq':-1, 's':-1, 'MIC':-1, 'ancestors':'-1'}
                return results

        offspring = pd.Series([birth_offspring(pop_sizes[k-1], w[k-1]) for k in range(1, n_species+1)], index=[f'{i}' for i in range(1, n_species+1)])

        mutants = pd.Series([get_mutants(x) for x in species], index=[i for i in species.index], dtype=object)

        for k in range(1, n_species+1):
            if len(mutants[f'{k}']) > 0:
                genotype_names = [i for i in species[f'{k}']['genotypes'].index]
                max_genotype_name = max([int(x.split('_')[-1]) for x in genotype_names])

                new_genotypes = [f'{current_prefix}_{max_genotype_name + i}' for i in range(1, len(mutants[f'{k}'])+1)]
                ancestor_genotype_names = pd.Series([species[f'{k}']['genotypes'][offspring[f'{k}'][i]]['name'] for i in mutants[f'{k}']])

                species[f'{k}'] = make_mutants(species[f'{k}'], ancestor_genotype_names, new_genotypes, mutant_func)
                offspring[f'{k}'][mutants[f'{k}']] = [len(genotype_names)-1+i for i in range(1, len(new_genotypes)+1)]
                
        genotype_ns = pd.Series([[sum(offspring[f'{k}'] == x) for x in range(len(species[f'{k}']['genotypes']))] for k in range(1, n_species+1)])
        species = pd.Series([set_n(species[f'{k}'], genotype_ns[k-1]) for k in range(1, n_species+1)], index=[f'{k}' for k in range(1, n_species+1)])
        for k in range(1, n_species+1):
            results[f'{k}'] = pd.concat((results[f'{k}'], pd.DataFrame({'gen':gen, 'genotypes':[i for i in species[f'{k}']['genotypes'].index], 'freq':get_n(species[f'{k}'])/species[f'{k}']['N'], 's':get_s(species[f'{k}']), 'MIC':get_MIC(species[f'{k}']), 'ancestors':get_ancestors(species[f'{k}'])})), ignore_index=True)

    return results

def summarize_results(results, wells, n_species, season):
    final = pd.DataFrame()
    for well in range(1, wells+1):
        for k in range(1, n_species+1):
            curr_result = results[f'{well}'][f'{k}']
            if any(curr_result['gen'] == -1):
                final = pd.concat((final, pd.DataFrame({'well':well, 'species':k, 'season':season, 'alive':False, 'n_genotypes':0, 's':np.nan, 'MIC':np.nan, 'mutations':np.nan}, index=[0])), ignore_index=True)
                continue
            curr_result = curr_result.loc[curr_result['gen']==max(curr_result['gen'])].loc[curr_result['freq']>0]
            curr_result = curr_result.assign(n_mutations=[(len(x.split(' '))-1) for x in curr_result['ancestors']])

            final = pd.concat((final, pd.DataFrame({'well':well, 'species':k, 'season':season, 'alive':True, 'n_genotypes':len(curr_result['genotypes']), 's':sum(curr_result['s']*curr_result['freq']), 'MIC':sum(curr_result['MIC']*curr_result['freq']), 'mutations':sum(curr_result['n_mutations']*curr_result['freq'])}, index=[0])), ignore_index=True)
    return final                

##### adamowicz_et_al_evolution_model_code_example.r
### simulation parameters
reps = 2 # reps per treatment condition ##
N = 100 # individuals per species ##
u = 0.01 # mutation rate ##
max_interdependent_species = 3 #
seasons = 5 # number of transfers ##
gens = 5 # gens per well per season ##
wells = 15 #
antibiotic_change_per_well = 1 # antibiotic concentration increase per well

### mutation function parameters
chance_of_MIC_mutation = 0.5 # probability that a mutation is a mutation in MIC
max_MIC_change = antibiotic_change_per_well * 1.1 # max mutation-induced MIC change
sd_s_change = 0.1 # if a mutation affects growth rate, the change is a normally-distributed variable centered on zero with this sd

### make mutation function
# null_function grabs MIC mutations from a uniform distribution
null_function = make_null_function(chance_of_MIC_mutation = chance_of_MIC_mutation, max_MIC_change = max_MIC_change, sd_s_change = sd_s_change)

### choose which mutation function to use
mutant_func = null_function
mutant_func_name = 'null'

### model
n_species_in_consortium = [i for i in range(1, max_interdependent_species+1)]
antibiotic = [i*antibiotic_change_per_well for i in range(wells)]
current_prefixes = [f'w{well}g1' for well in range(1, wells+1)]

all_data = pd.DataFrame()
for rep in range(1, reps+1):
    for n_species in n_species_in_consortium:
        '''
        print('rep', rep)
        print('\tgens', gens)
        print('\tu', u)
        print('\tmutant func', mutant_func_name)
        print('\tn_species', n_species)
        '''

        next_season = make_first_season(wells, n_species, current_prefixes, u = u)

        final = pd.DataFrame() ##or make it a list
        for season in range(seasons):
            results = [run_one_simulation(next_season[f'{well}'], gens, antibiotic[well-1], mutant_func, True, current_prefixes[well-1]) for well in [i for i in range(1, wells+1)]] ##
            next_season = start_next_season(results.copy(deep=True), wells, n_species, N, u=u)
            current_prefixes = [f'w{well}g{season+1}' for well in range(1, wells+1)]
            final = pd.concat((final, summarize_results(results, wells, n_species, season)), ignore_index=True)

        final = final.assign(n_species = n_species)
        final = final.assign(u = u)
        final = final.assign(rep = rep)
        final = final.assign(gens = gens)
        final = final.assign(mutant_function = mutant_func_name)
        all_data = pd.concat((all_data, final), ignore_index=True)
all_data.to_csv('all_data.csv', index=False)

###
all_data = pd.read_csv('all_data.csv', na_filter=False)

tol_data = all_data[['well', 'rep', 'gens', 'u', 'n_species', 'season', 'mutant_function']].loc[all_data['species'] == 1].loc[all_data['alive']].copy(deep=True)
#tol_data['n_species'] = tol_data['n_species'].astype(object)

tol_data = tol_data.groupby(['rep', 'gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).last()
tol_data = tol_data.assign(tolerance = tol_data['well'] - 1)
tol_stats = tol_data.groupby(['gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).tolerance.agg(['mean', 'std', 'count'])
tol_data = tol_data.groupby(['gens', 'u', 'n_species', 'season', 'mutant_function'], as_index=False).count()[['gens', 'u', 'n_species', 'season', 'mutant_function']]
tol_data = tol_data.assign(tolerance_sd = tol_stats['std'].reset_index(drop=True).copy(deep=True), n = tol_stats['count'].reset_index(drop=True).copy(deep=True), tolerance = tol_stats['mean'].reset_index(drop=True).copy(deep=True))

#sns.lineplot(x='season', y='tolerance', hue = 'n_species', err_style='bars', ci='sd', marker='o', data=tol_data)

error = tol_data['tolerance_sd'] / [math.sqrt(i-1) for i in tol_data['n']]

plt.errorbar(tol_data['season'].loc[tol_data['n_species']==1], tol_data['tolerance'].loc[tol_data['n_species']==1], error[:5], label = '1', color = 'tab:blue')#
plt.errorbar(tol_data['season'].loc[tol_data['n_species']==2], tol_data['tolerance'].loc[tol_data['n_species']==2], error[5:10], label = '2', color = 'tab:orange')#
plt.errorbar(tol_data['season'].loc[tol_data['n_species']==3], tol_data['tolerance'].loc[tol_data['n_species']==3], error[10:], label = '3', color = 'tab:green')#

plt.xlabel('transfer')
plt.ylabel('tolerance (arbitrary)')
plt.title(f'seed: {seed}')
plt.show()

## error bars, only 0.5's have tolerance_sd, tuples rather than dataframes
