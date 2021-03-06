# from https://youtu.be/MXUMJMrX2Gw

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t):
    '''
    x:list a list of outputs of each non-derived function at
    t:float time
    '''

    # constants
    a1 = 3e5
    a2 = 0.2
    a3 = 4e-7
    a4 = 0.6
    a5 = 8
    a6 = 90

    # assign each ODE to a vector element
    A = x[0]
    B = x[1]
    C = x[2]

    # define each ODE
    dAdt = a1 - a2*A - a3*A*C
    dBdt = a3*A*C - a4*B
    dCdt = -a3*A*C - a5*C + a6*B

    return [dAdt, dBdt, dCdt]

# define initial conditions
x0 = [2e6, 0, 90]

# (opt) test the defined ODEs
print(odes(x=x0, t=0))

# declare a time window
t = np.linspace(0, 15, 1000)
    ## endpoint=True

x = odeint(odes, x0, t)

A = x[:, 0]
B = x[:, 1]
C = x[:, 2]

# plot
plt.semilogy(t, A)
plt.semilogy(t, B)
plt.semilogy(t, C)
plt.show()