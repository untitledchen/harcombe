import math
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    rounded = math.floor(n * multiplier + 0.5) / multiplier
    if decimals == 0:
        rounded = int(rounded)
    return rounded

data = pd.read_csv('final_mono_132.csv', na_filter=False)
data = data.assign(name=data['genotype'] +  ", " +[str(round_half_up(i, 3)) for i in data['lag']])
data = data.assign(tp= [f'{i}, {j}' for i, j in zip(data['gen'], data['phase_end'])])
sns.relplot(data=data, kind='line', x='tp', y='ntot', hue='name')

sns.relplot(data=data.loc[data['phase_end'] == 2], kind='line', x='gen', y='ntot', hue='name')

data = pd.read_csv('final_mono_5.csv', na_filter=False)