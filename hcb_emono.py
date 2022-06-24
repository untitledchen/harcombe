from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import random
import math

### define odes equation
def odes(x, t, alpha, nE, tau_lag):
    '''
    x: list of initial conditions
    t: list of times
    alpha: additional death rate due to antibiotic, as a proportion of max growth rate r
    nE: number of strains of E. coli
    tau_lag: list of lag times of each strain # all identical for now
    '''
    
    M = x[0]
    L = x[1]

    # half-saturation constants
    K_M = 1 #
    K_L = 1 #

    # resource decay constants
    kM = 5e-9 #
    kL = 5e-9 #

    # E
    for i in range(1, nE + 1):

        # growth constants
        locals()[f'rE{i}'] = 1 #
        locals()[f'kE{i}'] = 5e-9 #
        locals()[f'tau_lagE{i}'] = tau_lag[i - 1]
        
        # resource constants
        locals()[f'cM{i}'] = 0.1 #
        locals()[f'cL{i}'] = 1.0 #

        # solutions -- starting with Eg1 on x[3] and El1 on x[4], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[2*i]
        locals()[f'El{i}'] = x[1 + 2*i]

        # differential equations
        locals()[f'dEg{i}dt'] = (1 - alpha)*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kE{i}']*locals()[f'Eg{i}'] + locals()[f'El{i}']/locals()[f'tau_lagE{i}']
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/locals()[f'tau_lagE{i}']

    # M
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dMdt = -sigma_cM - kM*M
    
    # L
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dLdt = -sigma_cL - kL*L

    to_return = [dMdt,dLdt]
    for i in range(1, nE + 1): ## could probably use list comprehension instead of this for loop
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])

    return to_return
###

### phase 1: antibiotic
def run_phase1(odes, init_cond, t_interval, nE, tau_lag):
    print(f'Running phase 1: antibiotic, with {init_cond}')
    
    sol = odeint(odes, init_cond, t_interval, args=(2, nE, tau_lag)) ## when alpha = -2, effective growth rate = -r
    para = [t_interval, nE, tau_lag]
    return sol, para

### phase 2: no antibiotic
def run_phase2(odes, sol1, t_interval, nE, tau_lag):
    
    init_cond = [sol1[:, 0][-1], sol1[:, 1][-1]]
    for i in range(1, nE + 1): ## could maybe use list comp
        init_cond += [sol1[:, (2*i)][-1], sol1[:, (1 + 2*i)][-1]]
    print(f'Running phase 2: no antibiotic, with {init_cond}')
    
    sol = odeint(odes, init_cond, t_interval, args=(0, nE, tau_lag))
    para = [t_interval, nE, tau_lag]
    return sol, para

### general matplotlib plot function
def solplot(sols, paras):
    '''
    sols:list of solutions for each phase
    t_intervals:list of t_intervals respective to sols
    '''
    fig, axs = plt.subplots(nrows=1, ncols=len(sols), figsize=(10, 4))

    for s in range(len(sols)):
        sol = sols[s]

        # read in parameters
        t_interval = paras[s][0]
        nE = paras[s][1]
        tau_lag = paras[s][2]
        
        # set left axis (nutrients) ylim to highest initial M or L overall (phase 1)
        nutr_ylim = max([max(sols[0][:, 0]), max(sols[0][:, 1])])
        
        locals()[f'M{s}'] = sol[:, 0]
        locals()[f'L{s}'] = sol[:, 1]
        for i in range(1, nE + 1):
            locals()[f'Eg{i}{s}'] = sol[:, (2*i)]
            locals()[f'El{i}{s}'] = sol[:, (1 + 2*i)]

        g = 0.8/(2*nE)
        for i in range(1, nE + 1):
            axs[s].semilogy(t_interval, locals()[f'Eg{i}{s}']+locals()[f'El{i}{s}'], color = (0, 0.8 - g*i, 0), label=f'E. coli strain {i}, lag time = {str(tau_lag[i - 1])}')

        axs[s].set_ylabel('log(Eg + El)')
        axs[s].set_xlabel('Time')
        
        locals()[f'axs{s}1'] = axs[s].twinx()
        
        locals()[f'axs{s}1'].plot(t_interval, locals()[f'M{s}'], label='methionine')
        locals()[f'axs{s}1'].plot(t_interval, locals()[f'L{s}'], label='lactose')
        
        locals()[f'axs{s}1'].set_ylim(0, nutr_ylim)

        locals()[f'axs{s}1'].set_ylabel('Amount of Nutrients')

        if s == 0: ## some hard code for titles and legends
            axs[s].set_title('Phase 1')
            axs[s].legend(loc='lower left')
            locals()[f'axs{s}1'].legend(loc='lower left', bbox_to_anchor=(0, 0.22))
        else:
            axs[s].set_title('Phase 2')

    fig.tight_layout()
    plt.show()

### run regular: observe as is
def regular(init_M, init_L, Ta, nE=1, tau_lag_inp=None):
    '''
    init_M: initial methionine concentration
    init_L: initial lactose concentration
    Ta: length of phase 1
    tau_lag: lag time, shared by all strains; or (default) lag time for each strain is a random value between 0 and 12 inclusive
    '''
    if tau_lag_inp == None:
        tau_lag = [random.uniform(0, 12) for i in range(nE)]
    else:
        tau_lag = [tau_lag_inp for i in range(nE)]
        
    init_cond = [init_M, init_L] + [0, 1]*nE ## each strain starts with growing population 0 and lagging population 1

    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, 15, 1000) #

    # run phases
    sol1, para1 = run_phase1(odes, init_cond, t_interval1, nE, tau_lag)
    
    sol2, para2 = run_phase2(odes, sol1, t_interval2, nE, tau_lag)

    #plot
    solplot([sol1, sol2], [para1, para2])
    
### run optimized: adapted to Fridman analysis
###
regular(10, 10, 6, nE=3)
