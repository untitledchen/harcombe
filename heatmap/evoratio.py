import pdb

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, linregress

def run_evoratio(filename_mono, filename_co, yvar):
    data1 = pd.read_csv(filename_mono, header=2, na_filter=False)
    data2 = pd.read_csv(filename_co, header=2, na_filter=False)

    data1[f'log10_{yvar}'] = np.log10(data1[yvar])
    data2[f'log10_{yvar}'] = np.log10(data2[yvar])

    slope_mono, intercept, r_value, p_value, std_err = linregress(data1['cycle'], data1[f'log10_{yvar}'])
    slope_co, intercept, r_value, p_value, std_err = linregress(data2['cycle'], data2[f'log10_{yvar}'])

    if slope_mono == 0:
        return -1
    return slope_co / slope_mono