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
        return [s, MIC]
    return null_function
            
##### model_functions2.r
def make_first_season(wells, n_species, current_prefixes, u=0.0025):
    next_season = pd.Series()
    for well in range(1, wells+1):
        current_prefix = current_prefixes[well]
        first_genotype = pd.Series([f'{current_prefix}_0', N, 1, 0, '0'], index=['name', 'n', 's', 'MIC', 'ancestors'])

        species = pd.Series()
        for k in range(1, n_species+1):
            species[k] = Species(N, u, first_genotype) ##

### simulation parameters
reps = 5 # reps per treatment condition
N = 1000 # individuals per species
u = 0.001 # mutation rate 
max_interdependent_species = 3 #
seasons = 20 # number of transfers
gens = 20 # gens per well per season
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
        print('rep', rep)
        print('\tgens', gens)
        print('\tu', u)
        print('\tmutant func', mutant_func_name)
        print('\tn_species', n_species)

        next_season = make_first_season(wells, n_species, current_prefixes, u = u)


        
