from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# for tweaking purposes
Ta = int(input('Ta: '))
tau_lag = int(input('tau_lag: '))
tlim2 = int(input('Time for phase 2: '))

def odes(x, t, a, tau_lag):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    # assign each ODE to a vector element
    G = x[0]
    L = x[1]

    # define each ODE
    dGdt = a*G + L/tau_lag
    dLdt = -L/tau_lag

    return [dGdt, dLdt]

# define initial conditions
x0 = [0, 1]

# (opt) test the defined ODEs
#print(odes(x=x0, t=0))

# declare a time window
t = np.linspace(0, Ta, 1000)
    # endpoint=True

x = odeint(odes, x0, t, args=(-1, tau_lag))

G = x[:, 0]
L = x[:, 1]

# plot
#plt.semilogy(t, G, label='growing')
#plt.semilogy(t, L, label='lagging')
#plt.legend()
plt.semilogy(t, G+L)
plt.title('antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.show()

###
x2 = [G[-1], L[-1]]
print(x2)
z = np.linspace(0, tlim2, 1000)

y = odeint(odes, x2, z, args=(1,tau_lag))

G2 = y[:, 0]
L2 = y[:, 1]

plt.semilogy(z, G2+L2)

plt.title('non-antibiotic phase')
plt.ylabel('log(g+l)')
plt.xlabel('t')
plt.ylim(sum(x2)*1e-2)

plt.show()

