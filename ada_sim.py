import pandas as pd
import random

##### mutation_functions.r
def make_null_function(chance_of_MIC_mutation, max_MIC_change, sd_s_change):
    def null_function():
        if random.uniform(0,1) < chance_of_MIC_mutation:
            MIC = random.uniform(0,1) * max_MIC_change # random proportion of max MIC change
            s = 0
        else:
            MIC = 0
            s = random.normalvariate(0, sd_s_change)
        return pd.Series([s, MIC])
    return null_function
            
##### model_functions2.r

# made class Species which is basically pd.Series but under a diff name
class SpeciesType(pd.core.series.Series):
    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False):
        pd.core.series.Series.__init__(self, data, index, dtype, name, copy, fastpath)

def Species(N=1000, u=0.001, first_genotype=None):
    me = pd.Series([pd.Series(dtype=object), N, u], index=['genotypes', 'N', 'u'], dtype=object)

    if type(first_genotype) != type(None): # avoiding ValueError "truth value of a Series is ambiguous"
        if {'name', 'n', 's', 'MIC', 'ancestors'}.issubset(set(first_genotype.index)):
            me['genotypes'][first_genotype['name']] = first_genotype
        else:
            print('All genotypes must contain the following fields: name (character), n (integer), s (number), MIC (number), ancestors (character)')
            return 0

    me = SpeciesType(me)
    return me

def get_s(Species):
    if type(Species) == SpeciesType:
        return pd.Series([x['s'] for x in Species['genotypes']])
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_n(Species):
    if type(Species) == SpeciesType:
        return pd.Series([x['n'] for x in Species['genotypes']])
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_MIC(Species):
    if type(Species) == SpeciesType:
        return pd.Series([x['MIC'] for x in Species['genotypes']])
    else:
        print('This method requires a SpeciesType object')
        return Species

def get_ancestors(Species):
    if type(Species) == SpeciesType:
        return pd.Series([x['ancestors'] for x in Species['genotypes']])
    else:
        print('This method requires a SpeciesType object')
        return Species

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
            print(sum(n))
            print('Sum of n != N defined in Species.')
            return -1
        for i in range(1, len(Species['genotypes'])):
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
            mutant = pd.Series({'name':new_genotype_names[k], 'n':1, 's':max([0, ancestor['s']+ mutant_function()[0]]), 'MIC':max([0, ancestor['s']+ mutant_function()[1]]), 'ancestors':ancestor_name+' '+ancestor['ancestors']})
            Species = add_genotype(Species, mutant)
            return Species
    else:
        print('This method requires a SpeciesType object')
        return Species

def make_first_season(wells, n_species, current_prefixes, u=0.0025):
    next_season = pd.Series(dtype=object)
    for well in range(1, wells+1):
        current_prefix = current_prefixes[well-1]
        first_genotype = pd.Series([f'{current_prefix}_0', N, 1, 0, '0'], index=['name', 'n', 's', 'MIC', 'ancestors'])
        
        species = pd.Series(dtype=object)
        for k in range(1, n_species+1):
            species[f'{k}'] = Species(N, u, first_genotype)

        next_season[f'{well}'] = species

    return next_season

def calculate_fitness(s, n, MIC, antibiotic):
    s[MIC < antibiotic] = 0
    
    if all(s==0):
        print('All individuals have zero fitness')
        return -1
    
    w = s/sum(s)*n
    w = w/sum(w)
    return w.cumsum()

def birth_offspring(N, w):
    chance = pd.Series([random.uniform(0,1) for i in range(N)])
    return pd.Series([w[x < w].index[0] for x in chance]) ## differs from code in that starts w 0

def run_one_simulation(species, gens, antibiotic, mutant_func, interdependent=True, current_prefix=''):
    n_species = len(species)
    pop_sizes = pd.Series([x['N'] for x in species])

    ## removed the len check
    antibiotic = pd.Series([antibiotic]*gens)

    if interdependent:
        if any(pd.Series([len(x['genotypes']) for x in species]) == 0):
            for k in range(1, n_species+1):
                results[f'{k}'] = pd.DataFrame({'gen':-1, 'genotypes':'-1', 'freq':-1, 's':-1, 'MIC':-1, 'ancestors':'-1'}, index=[0])
            return results

    results = pd.Series(dtype=object)
    for k in range(1, n_species+1):
        results[f'{k}'] = pd.DataFrame({'gen':0, 'genotypes': [i for i in species[f'{k}']['genotypes'].index], 'freq': get_n(species[f'{k}'])/species[f'{k}']['N'], 's': get_s(species[f'{k}']), 'MIC': get_MIC(species[f'{k}']), 'ancestors': get_ancestors(species[f'{k}'])}, index=[0])
    globals()['results'] = results#
    for gen in range(1, gens+1):
        w = pd.Series([calculate_fitness(get_s(k), get_n(k), get_MIC(k), antibiotic[gen-1]) for k in species])
        if interdependent:
            if any(w.explode() == -1): ##
                for k in range(1, n_species+1):
                    results[f'{k}'] = pd.concat((results[f'{k}'], pd.DataFrame({'gen':-1, 'genotypes':'-1', 'freq':-1, 's':-1, 'MIC':-1, 'ancestors':'-1'}, index=[0])), ignore_index=True)
                return results

        offspring = pd.Series([birth_offspring(pop_sizes[k-1], w[k-1]) for k in range(1, n_species+1)], index=[f'{i}' for i in range(1, n_species+1)])
        globals()['offspring'] = offspring#

        mutants = pd.Series([get_mutants(x) for x in species], index=[i for i in species.index], dtype=object)
        globals()['mutants'] = mutants#

        for k in range(1, n_species+1):
            if len(mutants[f'{k}']) > 0:
                genotype_names = [i for i in species[f'{k}']['genotypes'].index]
                max_genotype_name = max([int(x.split('_')[-1]) for x in genotype_names])

                new_genotypes = [f'{current_prefix}_{max_genotype_name + i}' for i in range(1, len(mutants[f'{k}'])+1)]

                ancestor_genotype_names = pd.Series([species[f'{k}']['genotypes'][offspring[f'{k}'][i]]['name'] for i in range(len(mutants[f'{k}']))])

                species[f'{k}'] = make_mutants(species[f'{k}'], ancestor_genotype_names, new_genotypes, mutant_func)
                offspring[f'{k}'][mutants[f'{k}']] = [len(genotype_names) + i for i in range(1, len(new_genotypes)+1)]

        genotype_ns = pd.Series([[sum(offspring[f'{k}'] == x) for x in range(len(species[f'{k}']['genotypes']))] for k in range(1, n_species+1)])
        species = pd.Series([set_n(species[f'{k}'], genotype_ns[k-1]) for k in range(1, n_species+1)])

        for k in range(1, n_species+1):
            print(k)
            results[f'{k}'] = pd.concat((results[f'{k}'], pd.DataFrame({'gen':gen, 'genotypes':[i for i in species[f'{k}']['genotypes'].index], 'freq':get_n(species[f'{k}'])/species[f'{k}']['N'], 's':get_s(species[f'{k}']), 'MIC':get_MIC(species[f'{k}']), 'ancestors':get_ancestors(species[f'{k}'])}, index=[0])), ignore_index=True)

    return results

### simulation parameters
reps = 1 # reps per treatment condition ##
N = 1000 # individuals per species
u = 0.001 # mutation rate 
max_interdependent_species = 2 # ##
seasons = 2 # number of transfers ##
gens = 3 # gens per well per season ##
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

numbers_of_species = n_species_in_consortium
antibiotic = [i*antibiotic_change_per_well for i in range(wells)]
current_prefixes = [f'w{well}g1' for well in range(1, wells+1)]

all_data = pd.DataFrame()
for rep in range(1, reps+1):
    for n_species in numbers_of_species:
        '''
        print('rep', rep)
        print('\tgens', gens)
        print('\tu', u)
        print('\tmutant func', mutant_func_name)
        print('\tn_species', n_species)
        '''

        next_season = make_first_season(wells, n_species, current_prefixes, u = u)

        final = pd.DataFrame()
        for season in range(1, seasons+1):
            results = pd.Series([run_one_simulation(next_season[f'{x}'], gens, antibiotic[x-1], mutant_func, True, current_prefixes[x-1]) for x in [i for i in range(1, wells + 1)]])

