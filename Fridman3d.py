from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math

# arbitrary
Talim = 11 # max Tas to try
tlim2 = 100 # max calculated t for n*0
taulim = 11 # max tau_lag to try per Ta

opt_taus = []
Tas = [i for i in range(1, Talim + 1)]
for Ta in Tas:

    # ode function
    def odes(x, t, a, tau_lag):

        G = x[0]
        L = x[1]

        dGdt = a*G + L/tau_lag
        dLdt = -L/tau_lag

        return [dGdt, dLdt]

    # n*0 (approximation) over tau_lag
    n_0 = []
    taus = np.linspace(0, taulim, 1000)[1:]

    for tau in taus:
        #phase 1
        x0 = [0, 1]
        t = np.linspace(0, Ta, 1000)
        x = odeint(odes, x0, t, args=(-1, tau))

        G = x[:, 0]
        L = x[:, 1]

        #phase 2
        x2 = [G[-1], L[-1]]
        z = np.linspace(0, tlim2, 1000)
        y = odeint(odes, x2, z, args=(1, tau))

        G2 = y[:, 0]
        L2 = y[:, 1]

        total = G2 + L2

        #n_0
        n_0.append((math.e)**-tlim2 * (total[-1]))

    # opt tau
    opt_taus.append(taus[n_0.index(max(n_0))])

plt.plot(Tas, opt_taus)
plt.ylabel('optimal tau lag')
plt.xlabel('Ta')
plt.show()
