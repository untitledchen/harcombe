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
'''
plt.semilogy(t, G+L)
plt.title('antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.show()
'''
### phase 2: no antibiotic, a = 1
x2 = [G[-1], L[-1]]
print(x2)
z = np.linspace(0, tlim2, 1000)

y = odeint(odes, x2, z, args=(1, tau_lag))

G2 = y[:, 0]
L2 = y[:, 1]
'''
plt.semilogy(z, G2+L2)

plt.title('non-antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.ylim(sum(x2)*1e-2)

plt.show()

'''
'''
### see what the calculated values of n*0 look like
total2 = G2+L2
n_0 = []
for i in range(len(z)):
    n_0.append((math.e)**(-(z[i])) * (total2[i]))
    
plt.plot(z, n_0)
plt.show()
'''
'''
# as t --> infinity (& slope approaches growth rate), calculated n*0 levels out to a certain value below 0.15

dGdt = []
for i in range(len(z)):
    dGdt.append(1*G2[i] + L2[i]/tau_lag)
    
dLdt = []
for i in range(len(z)):
    dLdt.append(-L2[i]/tau_lag)

dtotal = []
for i in range(len(dLdt)):
    dtotal.append(dGdt[i] + dLdt[i])

plt.semilogy(z, G2+L2, label='g+l')
plt.semilogy(z, dtotal, label='dg+dt')
plt.legend()
plt.show()
'''

###
n_0 = []
taus = np.linspace(0, 12, 1000)[1:]

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

    '''
    n_0 = []
    for i in range(len(z)):
        n_0.append((math.e)**(-(z[i])) * (total[i]))
    
    plt.plot(z, n_0)
    plt.show()
    '''
    
plt.plot(taus, n_0)
plt.show()
