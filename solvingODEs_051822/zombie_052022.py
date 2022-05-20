from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    # constants
    b = 0.5 # birth rate
    d = 0.4 # death rate
    beta = 0.2 # between-host transmission rate
    a = 0.4 # disease-induced death rate
    gamma = 0.1 # loss of immunity rate
    v = 0.1 # recovery rate

    # assign each ODE to a vector element
    S = x[0] # susceptible population
    I = x[1] # infected population
    R = x[2] # resistant population

    # define each ODE
    dSdt = b*S + b*I + b*R - d*S - beta*S*I + gamma*R
    dIdt = beta*S*I - (d + a)*I - v*I
    dRdt = v*I - d*R - gamma*R

    return [dSdt, dIdt, dRdt]

# define initial conditions
x0 = [5, 10, 6]

# (opt) test the defined ODEs
print(odes(x=x0, t=0))

# declare a time window
t = np.linspace(0, 15, 1000)
    ## endpoint=True

x = odeint(odes, x0, t)

S = x[:, 0]
I = x[:, 1]
R = x[:, 2]

# plot
plt.semilogy(t, S, label='S')
plt.semilogy(t, I, label='I')
plt.semilogy(t, R, label='R')
plt.legend()
plt.show()