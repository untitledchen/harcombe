import pdb

from hcb_sim import run
from calc_tolerance import run_calc_tol
import pandas as pd
import seaborn as sns
import random

import pdb

# seed = random.randrange(1000)
# print("running")
# run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
# run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))

data = pd.read_csv("C:\\Users\\untit\\harcombe\\data_hold\\hcb_sim_co_425_met1_lac1000.csv", header=1)

pdb.set_trace()
data = data.loc[data["phase_end"] == 2]
data["ntot"] = data["ngrow"] + data["nlag"]
data = data[["rep", "cycle", "species", "ntot"]]

data = data.groupby(["rep", "cycle", "species"]).sum()

# sns.hisplot(data, x="cycle", hue="species", multiple="stack")
# plt.show()

# data = pd.melt(data, id_vars=["rep", "cycle"], value_vars=["M", "L", "A"])

