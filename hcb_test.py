from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import math
import random

def odes(x, t, antibiotic_yn, nE=1, nS=1):
    '''
    antibiotic_yn: y = 1, n = -1
    '''
    
    M = x[0]
    A = x[1]
    L = x[2]

    # half-velocity constants
    K_M = 1 #
    K_L = 1 #
    K_A = 1 #

    # resource decay constants
    kM = 5e-9 #
    kA = 5e-9 #
    kL = 5e-9 #

    # E
    for i in range(1, nE + 1):

        # constants
        locals()[f'rE{i}'] = 1*antibiotic_yn #
        locals()[f'kE{i}'] = 5e-9 #
        globals()[f'tau_lagE{i}'] =  random.uniform(1,12) #

        # resource constants
        locals()[f'cM{i}'] = 0.1 #
        locals()[f'pA{i}'] = 1.01 #
        locals()[f'cL{i}'] = 1.0 #

        # solutions -- starting with Eg1 on x[3] and El1 on x[4], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[1 + 2*i]
        locals()[f'El{i}'] = x[2 + 2*i]

        # differential equations
        locals()[f'dEg{i}dt'] = locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kE{i}']*locals()[f'Eg{i}'] + locals()[f'El{i}']/globals()[f'tau_lagE{i}'] #
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/globals()[f'tau_lagE{i}'] #

    # S
    for j in range(1, nS + 1):

        # constants
        locals()[f'rS{j}'] = 0.5*antibiotic_yn #
        locals()[f'kS{j}'] = 5e-9 #
        locals()[f'tau_lagS{j}'] = 6 # fixed

        # resource constants
        locals()[f'cA{j}'] = 1.0 #
        locals()[f'pM{j}'] = 1.56 #
        
        # solutions -- starting after the last strain of El, Sg and Sl occupy alternating indexes for each strain
        locals()[f'Sg{j}'] = x[(1 + 2*nE) + 2*j]
        locals()[f'Sl{j}'] = x[(2 + 2*nE) + 2*j]

        # differential equations
        locals()[f'dSg{j}dt'] = locals()[f'rS{j}']*locals()[f'Sg{j}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kS{j}']*locals()[f'Sg{j}'] + locals()[f'Sl{j}']/locals()[f'tau_lagS{j}']
        locals()[f'dSl{j}dt'] = -locals()[f'Sl{j}']/locals()[f'tau_lagS{j}']

    # M
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}']*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
    sigma_pM = 0
    for j in range(1, nS + 1):
        sigma_pM += locals()[f'pM{j}']*locals()[f'rS{j}']*locals()[f'Sg{j}']*(A/(A+K_A))
        
    dMdt = -sigma_cM + sigma_pM - kM*M
    
    # A
    sigma_pA = 0
    for i in range(1, nE + 1):
        sigma_pA += locals()[f'pA{i}']*locals()[f'rE{i}']*locals()[f'Eg{i}']*((M+0.1*K_M)/(M+K_M))*(L/(L+K_L)) #note kick-off
    sigma_cA = 0
    for j in range(1, nS + 1):
        sigma_cA += locals()[f'cA{j}']*locals()[f'rS{j}']*locals()[f'Sg{j}']*(A/(A+K_A))
        
    dAdt = sigma_pA - sigma_cA - kA*A
    
    # L
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}']*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dLdt = -sigma_cL - kL*L

    to_return = [dMdt, dAdt, dLdt]
    for i in range(1, nE + 1): # could probably use list comprehension instead of this for loop
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])
    for j in range(1, nS + 1):
        to_return.append(locals()[f'dSg{j}dt'])
        to_return.append(locals()[f'dSl{j}dt'])

    return to_return

def run_phase1(init_cond, odes, t_interval, Ta, nE, nS):

    sol1 = odeint(odes, init_cond, t_interval, args=(1, nE, nS))

    return sol1

def run_phase2(sol):
    end_cond = []
    for n in len(sol): # again, there's probably a better way than a loop
        end_cond.append(sol[:,n])

    tlim = 100 #
    t_interval2 = np.linspace(0, tlim, 1000)

    sol2 = odeint(odes, end_cond, t_interval2, args=(-1, nE, nS))

    return sol2

def start():
    init_cond = [1, 1, 1, 0, 1, 0, 1] #
    Ta = 6 #
    nE = 1 #
    nS = 1 #

    # phase 1
    print('Running phase 1: antibiotic present. Growth rate = maximum growth rate * -1')
    print('Initial conditions:')
    print(f'Ta: {Ta}')
    print(f'M: {init_cond[0]}')
    print(f'A: {init_cond[1]}')
    print(f'L: {init_cond[2]}')
    for i in range(1, nE + 1):
        print(f'Eg{i}: {init_cond[1 + 2*i]}')
        print(f'El{i}: {init_cond[2 + 2*i]}')
    for j in range(1, nS + 1):
        print(f'Sg{j}: {init_cond[(1 + 2*nE) + 2*j]}')
        print(f'Sl{j}: {init_cond[(2 + 2*nE) + 2*j]}')

    t_interval1 = np.linspace(0, Ta, 1000)
    sol1 = run_phase1(init_cond, odes, t_interval1, Ta, nE, nS)
    M = sol1[:, 0]
    A = sol1[:, 1]
    L = sol1[:, 2]
    for i in range(1, nE + 1):
        locals()[f'Eg{i}'] = sol1[:, (1 + 2*i)]
        locals()[f'El{i}'] = sol1[:, (2 + 2*i)]
    for j in range(1, nS + 1):
        locals()[f'Sg{j}'] = sol1[:, (1 + 2*nE) + 2*j]
        locals()[f'Sl{j}'] = sol1[:, (2 + 2*nE) + 2*j]

    for i in range(1, nE + 1):
        plt.semilogy(t_interval1, locals()[f'Eg{i}']+locals()[f'El{i}'], label=str(globals()[f'tau_lagE{i}'])) #
    for j in range(1, nS + 1):
        plt.semilogy(t_interval1, locals()[f'Sg{j}']+locals()[f'Sl{j}'])
    plt.ylabel('log(Eg + El)')
    plt.xlabel('t')
    plt.legend()

    plt.show()
    
