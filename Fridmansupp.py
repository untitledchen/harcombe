from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math

# for tweaking purposes
Ta = 6
tau_lag = 6
tlim2 = 15

### phase 1: antibiotic, a = -1
def odes(x, t, a, tau_lag):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    G = x[0]
    L = x[1]

    dGdt = a*G + L/tau_lag
    dLdt = -L/tau_lag

    return [dGdt, dLdt]

# initial conditions
x0 = [0, 1]

t = np.linspace(0, Ta, 1000)

x = odeint(odes, x0, t, args=(-1, tau_lag))

G = x[:, 0]
L = x[:, 1]

plt.semilogy(t, G+L)
plt.title('antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.show()

### phase 2: no antibiotic, a = 1
x2 = [G[-1], L[-1]]
print('[G, L] at end of phase 1:', x2)
z = np.linspace(0, tlim2, 1000)

y = odeint(odes, x2, z, args=(1, tau_lag))

G2 = y[:, 0]
L2 = y[:, 1]

plt.semilogy(z, G2+L2)

plt.title('non-antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.ylim(sum(x2)*1e-2)

plt.show()

### n*0 (approximation) over tau_lag
n_0 = []
taus = np.linspace(0, 12, 1000)[1:] # to avoid divide by zero error

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

plt.plot(taus, n_0)
plt.ylabel('n*0')
plt.xlabel('tau lag')
plt.show()

# testing to find max
print(max(n_0))
print(n_0.index(max(n_0)))
print(n_0[n_0.index(max(n_0))])
print(taus[n_0.index(max(n_0))])
