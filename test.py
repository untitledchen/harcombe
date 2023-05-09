import numpy as np
import time

for i in range(3):
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    z = np.array(y)
    begin = time.perf_counter()
    b = list(x).extend(y)
    print(time.perf_counter() - begin)
    a = np.concatenate((x, z), axis=0)
    print(time.perf_counter() - begin)

    # 1.09999964479357e-05
    # 7.270000060088933e-05
    # 6.300004315562546e-06
    # 2.649999805726111e-05
    # 4.899993655271828e-06
    # 2.2599997464567423e-05