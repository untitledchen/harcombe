# from https://youtu.be/zRMmiBMjP9o

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    # constants
    kr1 = 1e5 # new healthy cells a year
    kr2 = 0.1 # death rate of healthy cells
    kr3 = 2e-7 # healthy cell conversion rate to infected
    kr4 = 0.5 # death rate of infected cells
    kr5 = 5 # death rate of virus
    kr6 = 100 # production rate of virus by infected cells

    # assign each ODE to a vector element
    H = x[0] # healthy cells
    I = x[1] # infected cells
    V = x[2] # virus count

    # define each ODE
    dHdt = kr1 - kr2*H - kr3*H*V
    dIdt = kr3*H*V - kr4*I
    dVdt = -kr3*H*V - kr5*V + kr6*I

    return [dHdt, dIdt, dVdt]

# define initial conditions
x0 = [1e6, 0, 100]

# (opt) test the defined ODEs
print(odes(x=x0, t=0))

# declare a time window
t = np.linspace(0, 15, 1000)
    ## endpoint=True

x = odeint(odes, x0, t)

H = x[:, 0]
I = x[:, 1]
V = x[:, 2]

# plot
plt.semilogy(t, H, label='healthy cells')
plt.semilogy(t, I, label='infected cells')
plt.semilogy(t, V, label='virus count')
plt.legend()
plt.show()