from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math
import copy

### define odes function
def odes(x, t, alpha, nE, tau_lag):
    '''
    alpha: additional death rate due to antibiotic, as a proportion of max growth rate r
    '''
    
    M = x[0]
    L = x[1]

    # half-velocity constants
    K_M = 1 #
    K_L = 1 #
    K_A = 1 #

    # resource decay constants
    kM = 5e-9 #
    kL = 5e-9 #

    # E
    for i in range(1, nE + 1):

        # constants
        locals()[f'rE{i}'] = 1 #
        locals()[f'kE{i}'] = 5e-9 #
        globals()[f'tau_lagE{i}'] = tau_lag ## set all strains to same tau_lag for now

        # resource constants
        locals()[f'cM{i}'] = 0.1 #
        locals()[f'cL{i}'] = 1.0 #

        # solutions -- starting with Eg1 on x[3] and El1 on x[4], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[2*i]
        locals()[f'El{i}'] = x[1 + 2*i]

        # differential equations
        locals()[f'dEg{i}dt'] = (1 - alpha)*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kE{i}']*locals()[f'Eg{i}'] + locals()[f'El{i}']/globals()[f'tau_lagE{i}'] ## note globals() usage
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/globals()[f'tau_lagE{i}'] ## note globals() usage

    # M
    '''
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
    ''' 
    dMdt = 10
    
    # L
    '''
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
    '''  
    dLdt = 10

    to_return = [dMdt,dLdt]
    for i in range(1, nE + 1): ## could probably use list comprehension instead of this for loop
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])

    return to_return
###

def run_phase1(odes, init_cond, t_interval, nE, tau_lag):

    sol1 = odeint(odes, init_cond, t_interval, args=(2, nE, tau_lag)) ## when alpha = -2, effective growth rate = -r
    return sol1

def run_phase2(odes, init_cond, t_interval, nE, tau_lag):

    sol2 = odeint(odes, init_cond, t_interval, args=(0, nE, tau_lag))
    return sol2

### n*0 (approximation) over tau_lag
tlim2 = 15 #
t_interval2 = np.linspace(0, tlim2, 1000)

nE = 3 #
init_cond = [10, 10] + [0, 1]*nE

Talim = 11 # max Tas to try
opt_taus = [] # list of tau_lags that produce the maximum n*0 for their specific Ta
Tas = [i for i in range(1, Talim + 1)] # list of Tas to try, whole numbers only
#print('M and L kept constant at 10')
for Ta in Tas:

    n_0 = [] # list of all calculated n*0s for each tau_lag tested
    taus = np.linspace(0, 12, 1000)[1:] # list of tau_lags to be tested, ranging from just above 0 to 12.00
    # [1:] to avoid a divide by zero error
    
    for tau in taus:

        # phase 1, antibiotics
        t_interval1 = np.linspace(0, Ta, 1000)
        sol1 = run_phase1(odes, init_cond, t_interval1, nE, tau)
            
        M = sol1[:, 0]
        L = sol1[:, 1]
        for i in range(1, nE + 1):
            locals()[f'Eg{i}'] = sol1[:, (2*i)]
            locals()[f'El{i}'] = sol1[:, (1 + 2*i)]

        # phase 2, no antibiotics
        init_cond2 = [M[-1], L[-1]]
        for i in range(1, nE + 1):
            init_cond2 += [locals()[f'Eg{i}'][-1], locals()[f'El{i}'][-1]]
            
        sol2 = run_phase2(odes, init_cond2, t_interval2, nE, tau)
            
        M = sol2[:, 0]
        L = sol2[:, 1]
        for i in range(1, nE + 1):
            locals()[f'Eg{i}'] = sol2[:, (2*i)]
            locals()[f'El{i}'] = sol2[:, (1 + 2*i)]

        # let's just focus on strain 1 for now
        total = Eg1 + El1

        # append calculated n*0 to n_0
        n_0.append((math.e)**-tlim2 * (total[-1])) # calculate n*0 using the last t value available (== tlim2)

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
