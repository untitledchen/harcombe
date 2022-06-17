from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math
import random

###
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

def run_phase1(odes, init_cond, t_interval, nE, tau_lag):

    sol1 = odeint(odes, init_cond, t_interval, args=(2, nE, tau_lag)) ## when alpha = -2, effective growth rate = -r
    return sol1

def run_phase2(odes, init_cond, t_interval, nE, tau_lag):

    sol2 = odeint(odes, init_cond, t_interval, args=(0, nE, tau_lag))
    return sol2

def start(init_M, init_L):
    print(f'Running start({init_M}, {init_L}). Initial methionine = {init_M}, initial lactose = {init_L}.')
    tau_lag = 6 #
    
    Ta = 6 #
    nE = 3 #
    init_cond = [init_M, init_L] + [0, 1]*nE

    # phase 1
    print('Running phase 1: antibiotic present. Growth rate = maximum growth rate * -1')
    print('Initial conditions:')
    print(f'Ta: {Ta}')
    print(f'M: {init_cond[0]}')
    print(f'L: {init_cond[1]}')
    for i in range(1, nE + 1):
        print(f'Eg{i}: {init_cond[2*i]}')
        print(f'El{i}: {init_cond[1 + 2*i]}')

    t_interval1 = np.linspace(0, Ta, 1000)
    sol1 = run_phase1(odes, init_cond, t_interval1, nE, tau_lag)
    
    M = sol1[:, 0]
    L = sol1[:, 1]
    for i in range(1, nE + 1):
        locals()[f'Eg{i}'] = sol1[:, (2*i)]
        locals()[f'El{i}'] = sol1[:, (1 + 2*i)]

    # plot
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    g = 0.8/(2*nE)
    for i in range(1, nE + 1):
        axs[0].semilogy(t_interval1, locals()[f'Eg{i}']+locals()[f'El{i}'], color = (0, 0.8 - g*i, 0), label=f"E. coli strain {i}, lag time = {str(globals()[f'tau_lagE{i}'])}") ## note globals() usage

    axs[0].set_title('Phase 1')
    axs[0].set_ylabel('log(Eg + El)')
    axs[0].set_xlabel('Time')
    axs[0].legend(loc='lower left')
    
    axs01 = axs[0].twinx()
    
    axs01.plot(t_interval1, M, label='methionine')
    axs01.plot(t_interval1, L, label='lactose')

    
    axs01.set_ylim(0, max([init_M, init_L]))

    axs01.set_ylabel('Amount of Nutrients')
    axs01.legend(loc='lower left', bbox_to_anchor=(0, 0.22))

    tlim2 = 15 #
    init_cond2 = [M[-1], L[-1]]
    for i in range(1, nE + 1): ## list comp?
        init_cond2 += [locals()[f'Eg{i}'][-1], locals()[f'El{i}'][-1]]
        
    # phase 2
    print('Running phase 2: antibiotic absent. Growth rate = maximum growth rate * 1')
    print('Initial conditions (phase 1 final conditions):')
    print(f'M: {M[-1]}')
    print(f'L: {L[-1]}')
    for i in range(1, nE + 1):
        print(f"Eg{i}: {locals()[f'Eg{i}'][-1]}")
        print(f"El{i}: {locals()[f'El{i}'][-1]}")

    t_interval2 = np.linspace(0, tlim2, 1000)
    sol2 = run_phase2(odes, init_cond2, t_interval2, nE, tau_lag)
    
    M = sol2[:, 0]
    L = sol2[:, 1]
    for i in range(1, nE + 1):
        locals()[f'Eg{i}'] = sol2[:, (2*i)]
        locals()[f'El{i}'] = sol2[:, (1 + 2*i)]

    # plot
    g = 0.8/(2*nE)
    for i in range(1, nE + 1):
        axs[1].semilogy(t_interval2, locals()[f'Eg{i}']+locals()[f'El{i}'], color = (0, 0.8 - g*i, 0), label=f"E. coli strain {i}, lag time = {str(globals()[f'tau_lagE{i}'])}") ## note globals() usage

    axs[1].set_title('Phase 2')
    axs[1].set_ylabel('log(Eg + El)')
    axs[1].set_xlabel('Time')
    
    axs11 = axs[1].twinx()
    
    axs11.plot(t_interval2, M, label='methionine')
    axs11.plot(t_interval2, L, label='lactose')
    
    axs11.set_ylabel('Amount of Nutrients')

    axs11.set_ylim(0, max([init_M, init_L]))

    fig.tight_layout()
    plt.show()
    
start(1, 1)

