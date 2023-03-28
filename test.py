import numpy as np

x = [0, 4, 2, -6, 4, -2, 10]
y = [[0, i][i > 0] for i in x]