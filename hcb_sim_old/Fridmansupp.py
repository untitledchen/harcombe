# replicate the figures shown in the Fridman supplementary materials

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math

### arbitrary constants
Ta = 6 # the duration of phase 1 (hr)
tau_lag = 6 # the lag time (hr)
tlim2 = 15 # effectively, the duration of phase 2. This max t value will be used to calculate n*0 (greater tlim2 -> greater accuracy of n*0) 

### define the ODE function
def odes(x, t, a, tau_lag):
    '''
    x: a list of outputs of each non-derived function at
    t: time
    a: growth constant
    tau_lag: lag time (hr)
    '''

    G = x[0]
    L = x[1]

    dGdt = a*G + L/tau_lag
    dLdt = -L/tau_lag

    return [dGdt, dLdt]

### phase 1: antibiotic (growth constant, a = -1)

x0 = [0, 1] # initial conditions
t_interval = np.linspace(0, Ta, 1000) # time interval

x = odeint(odes, x0, t_interval, args=(-1, tau_lag)) # solve ODEs
# and assign the solved values to variables G and L
G = x[:, 0]
L = x[:, 1]

# plot
plt.semilogy(t_interval, G+L)

plt.title('Antibiotic Phase, a = -1')
plt.ylabel('log(g+l)')
plt.xlabel('t')

plt.show()

### phase 2: no antibiotic, a = 1
x2 = [G[-1], L[-1]] # initial conditions: the final values of G and L from phase 1
print('[G, L] at end of phase 1:', x2)
t_interval2 = np.linspace(0, tlim2, 1000) # time interval

y = odeint(odes, x2, t_interval2, args=(1, tau_lag)) # solve ODEs
# and assign the solved values to variables G2 and L2
G2 = y[:, 0]
L2 = y[:, 1]

# plot
plt.semilogy(t_interval2, G2+L2)

plt.title('Non-antibiotic Phase, a = 1')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.ylim(sum(x2)*1e-2) # set the y-axis range to a reasonable distance below 0, so we can see where the asymptotic line would approximately intersect the y-axis

plt.show()

### n*0 (approximation) over tau_lag
n_0 = [] # list of all calculated n*0s for each tau_lag tested
taus = np.linspace(0, 12, 1000)[1:] # list of tau_lags to be tested, ranging from just above 0 to 12.00
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
    n_0.append((math.e)**-tlim2 * (total[-1])) # calculate n*0 using the last t value available (== tlim2)

# plot: tau_lag vs. calculated n*0
plt.plot(taus, n_0)

plt.title('Effective Population (n*0) as a Result of Lag Time (tau_lag)')
plt.ylabel('n*0')
plt.xlabel('tau_lag')

plt.show()

# testing to find max n*0 / optimal tau_lag
print('')
print('max n*0:', max(n_0))
print('index of n_0 in which max n*0 is found:', n_0.index(max(n_0)))
print(f'(double-checking to see that indexing n_0 with the previous value does yield the max n*0: {n_0[n_0.index(max(n_0))]})')
print('tau_lag value for max n*0, using the previous index on the list of tau_lag values:', taus[n_0.index(max(n_0))])
