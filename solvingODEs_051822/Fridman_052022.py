from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    # constants
    a = -1
    tau_lag = 1 ##

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
print(odes(x=x0, t=0))

# declare a time window
t = np.linspace(0, 15, 1000)
    ## endpoint=True

x = odeint(odes, x0, t)

G = x[:, 0]
L = x[:, 1]

# plot
plt.semilogy(t, G, label='growing')
plt.semilogy(t, L, label='lagging')
plt.legend()
plt.show()