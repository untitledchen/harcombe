import numpy as np
import pandas as pd
import time

for i in range(3):
    data = pd.read_csv('phase0/hcb_sim_co_166_met1_phase00.csv', header=1, na_filter=False)

    begin = time.perf_counter()
    b = len(data.loc[data['species'] == 'Escherichia coli'].index)
    print(time.perf_counter() - begin)
    a = len(data.loc[data['species'] == 'Escherichia coli'])
    print(time.perf_counter() - begin)

    # 1.09999964479357e-05
    # 7.270000060088933e-05
    # 6.300004315562546e-06
    # 2.649999805726111e-05
    # 4.899993655271828e-06
    # 2.2599997464567423e-05


    # data = pd.read_csv('hcb_sim_co_166_met1_phase00.csv', header=1, na_filter=False)
    #
    # begin = time.perf_counter()
    # b = len(data.loc[data['species'] == 'Escherichia coli'].index)
    # print(time.perf_counter() - begin)
    # a = len(data.loc[data['species'] == 'Escherichia coli'])
    #
    # 0.0016824000049382448
    # 0.0025371999945491552
    # 0.0011260999599471688
    # 0.002015499980188906
    # 0.0012128999806009233
    # 0.0020409999997355044
