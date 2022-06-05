# recreate Fig 3d from Fridman paper

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math

# arbitrary constants
Talim = 11 # max Tas to try
tlim2 = 100 # max calculated t for phase 2, used to find n*0
taulim = 11 # max tau_lag to try per Ta

opt_taus = [] # list of tau_lags that produce the maximum n*0 for their specific Ta
Tas = [i for i in range(1, Talim + 1)] # list of Tas to try, whole numbers only

for Ta in Tas:

    # ode function
    def odes(x, t, a, tau_lag):

        G = x[0]
        L = x[1]

        dGdt = a*G + L/tau_lag
        dLdt = -L/tau_lag

        return [dGdt, dLdt]

    # n*0 (approximation) over tau_lag
    n_0 = [] # list of all calculated n*0s for each tau_lag tested
    taus = np.linspace(0, taulim, 1000)[1:] # list of tau_lags to be tested, ranging from just above 0 to 12.00
        # [1:] to avoid a divide by zero error

    for tau in taus:
        # phase 1, antibiotics
        x0 = [0, 1]
        t_interval = np.linspace(0, Ta, 1000)
        x = odeint(odes, x0, t_interval, args=(-1, tau))

        G = x[:, 0]
        L = x[:, 1]

        # phase 2, no antibiotics
        x2 = [G[-1], L[-1]]
        t_interval2 = np.linspace(0, tlim2, 1000)
        y = odeint(odes, x2, t_interval2, args=(1, tau))

        G2 = y[:, 0]
        L2 = y[:, 1]

        total = G2 + L2

        # append calculated n*0 to n_0
        n_0.append((math.e)**-tlim2 * (total[-1])) # calculated n*0 using the last t value (== tlim2)

    # append the tau_lag corresponding to max(n_0) to opt_taus
    opt_taus.append(taus[n_0.index(max(n_0))])

# plot
fig = plt.figure()
ax = fig.add_subplot()

plt.plot(Tas, opt_taus)

plt.title('Optimal Lag Time (Optimal tau_lag) Corresponding to a Given Duration of Antibiotic Treatment (Ta)')
plt.ylabel('Optimal tau_lag')
plt.xlabel('Ta')
ax.set_aspect('equal', adjustable='box')

plt.show()
