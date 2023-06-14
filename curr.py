from hcb_sim import run
from calc_tolerance import run_calc_tol
import pandas as pd
import random

# seed = random.randrange(1000)
# print("running")
# run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
# run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
#
# run(seed, "mono", 10, (0.01, 0), 10, (80, 2780, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
# run(seed, "co", 10, (0.01, 0.01), 10, (1, 2780, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

run_calc_tol(input('input'), 1000, 0.01, 0.1, 0.5)
run_calc_tol(input('input'), 1000, 0.01, 0.1, 0.5)