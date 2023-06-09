from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import random
import math


##issue with odes
### define odes equation
def odes(x, t, alpha, tau_lag, n, frid=False):
    '''
    x:[M, A, L, Egi, Eli, Sgi, Sli] initial conditions where Eg and El alternate indexes for i strains, then are followed by Sg and Sl alternating for j strains
    t:array of times
    alpha:[alphaE, alphaS] where alphaE: additional death rate due to antibiotic, as a proportion of max growth rate r for E. coli, alphaS: "" for S.enterica
    tau_lag:[tau_lagE, tau_lagS] where tau_lagE:list of tau_lags for each strain of E. coli, tau_lagS: "" for S. enterica 
    
    n:[nE, nS] where nE:int strains of E. coli, nS: "" for S. enterica
    frid: when True, dMdt & dAdt & dLdt = 0 (for use with Fridman-style analysis)
    '''
    # read in parameters
    alphaE = alpha[0]
    alphaS = alpha[1]

    tau_lagE = tau_lag[0]
    tau_lagS = tau_lag[1]

    nE = n[0]
    nS = n[1]

    # R
    M = x[0]
    A = x[1]
    L = x[2]

    # half-saturation constants
    K_M = 1 #
    K_A = 1 #
    K_L = 1 #

    # resource decay constants
    kM = 5e-9 #
    kA = 5e-9 #
    kL = 5e-9 #

    # E
    for i in range(1, nE + 1):

        # growth constants
        locals()[f'rE{i}'] = 1 #
        locals()[f'kE{i}'] = 5e-9 #
        locals()[f'tau_lagE{i}'] = tau_lagE[i- 1]
        
        # resource constants
        locals()[f'cM{i}'] = 0.1 #
        locals()[f'pA{i}'] = 1.01 #
        locals()[f'cL{i}'] = 1.0 #

        # solutions -- starting with Eg1 on x[3] and El1 on x[4], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[1 + 2*i]
        locals()[f'El{i}'] = x[2 + 2*i]

        # differential equations
        locals()[f'dEg{i}dt'] = (1 - alphaE)*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kE{i}']*locals()[f'Eg{i}'] + locals()[f'El{i}']/locals()[f'tau_lagE{i}']
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/locals()[f'tau_lagE{i}']

    # S
    for j in range(1, nS + 1):

        # constants
        locals()[f'rS{j}'] = 0.5 #
        locals()[f'kS{j}'] = 5e-9 #
        locals()[f'tau_lagS{j}'] = tau_lagS[j - 1]

        # resource constants
        locals()[f'cA{j}'] = 1.0 #
        locals()[f'pM{j}'] = 1.56 #
        
        # solutions -- starting with Sg1 after the last El, Sg and Sl occupy alternating indexes for each strain
        locals()[f'Sg{j}'] = x[(1 + 2*nE) + 2*j]
        locals()[f'Sl{j}'] = x[(2 + 2*nE) + 2*j]

        # differential equations
        locals()[f'dSg{j}dt'] = (1 - alphaS)*locals()[f'rS{j}']*locals()[f'Sg{j}']*(A/(A+K_A)) - locals()[f'kS{j}']*locals()[f'Sg{j}'] + locals()[f'Sl{j}']/locals()[f'tau_lagS{j}']
        locals()[f'dSl{j}dt'] = -locals()[f'Sl{j}']/locals()[f'tau_lagS{j}']

    # M
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
    sigma_pM = 0
    for j in range(1, nS + 1):
        sigma_pM += locals()[f'pM{j}']*locals()[f'rS{j}']*locals()[f'Sg{j}']*(A/(A+K_A))
     
    dMdt = (-sigma_cM + sigma_pM - kM*M) * ([1, 0][frid])

    # A
    sigma_pA = 0
    for i in range(1, nE + 1):
        sigma_pA += locals()[f'pA{i}']*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
    sigma_cA = 0
    for j in range(1, nS + 1):
        sigma_cA += locals()[f'cA{j}']*locals()[f'Sg{j}']*(A/(A+K_A))
        
    dAdt = (sigma_pA - sigma_cA - kA*A) * ([1, 0][frid])
    
    # L
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dLdt = (-sigma_cL - kL*L) * ([1, 0][frid])

    to_return = [dMdt, dAdt, dLdt]
    for i in range(1, nE + 1):
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])
    for j in range(1, nS + 1):
        to_return.append(locals()[f'dSg{j}dt'])
        to_return.append(locals()[f'dSl{j}dt'])

    return to_return
###

### run_phase 1 (antibiotic) or 2 (no antibiotic)
def run_phase(odes, init_cond, t_interval, tau_lag, n, frid=False, phase=1):

    if phase == 1:
        alpha = [2, 2]
    elif phase == 2:
        alpha = [0, 0]
    
    sol = odeint(odes, init_cond, t_interval, args=(alpha, tau_lag, n, frid)) ## when alpha = -2, effective growth rate = -r
    para = [t_interval, tau_lag, n]
    return sol, para

### general matplotlib plot function
def solplot(sols, paras):
    '''
    sols:list of solutions for each phase
    paras:list of parameters for each phase [t_interval, tau_lag, n]
    '''
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))

    for s in range(len(sols)):
        sol = sols[s]

        # read in parameters
        t_interval = paras[s][0]
        tau_lag = paras[s][1]
        
        n = paras[s][2]
        nE = n[0]
        nS = n[1]
        
        locals()[f'M{s}'] = sol[:, 0]
        locals()[f'A{s}'] = sol[:, 1]
        locals()[f'L{s}'] = sol[:, 2]
        for i in range(1, nE + 1):
            locals()[f'Eg{i}{s}'] = sol[:, (1 + 2*i)]
            locals()[f'El{i}{s}'] = sol[:, (2 + 2*i)]
        for j in range(1, nS + 1):
            locals()[f'Sg{j}{s}'] = sol[:, ((1 + 2*nE) + 2*j)]
            locals()[f'Sl{j}{s}'] = sol[:, ((2 + 2*nE) + 2*j)]

        # create totals to convert absolute populations to frequencies of the total population, by creating totals and later freqs
        totals = locals()[f'Eg1{s}'] + locals()[f'El1{s}'] + locals()[f'Sg1{s}'] + locals()[f'Sl1{s}']
        for i in range(2, nE + 1):
            totals += locals()[f'Eg{i}{s}'] + locals()[f'El{i}{s}']
        for j in range(2, nS + 1):
            totals += locals()[f'Sg{j}{s}'] + locals()[f'Sl{j}{s}']

        g = 0.8/(1.5*nE)
        for i in range(1, nE + 1):
            freqsE = (locals()[f'Eg{i}{s}']+locals()[f'El{i}{s}'])/totals
            axs[0][s].plot(t_interval, np.log10(freqsE), color = (0, 0.8 - g*i, 0), label=f'E. coli strain {i}, lag time = {str(tau_lag[0][i - 1])}')
        b = 0.8/(1.5*nS)
        for j in range(1, nS + 1):
            freqsS = (locals()[f'Sg{j}{s}']+locals()[f'Sl{j}{s}'])/totals
            axs[0][s].plot(t_interval, np.log10(freqsS), color = (0, 0, 0.8 - b*j), label=f'S. enterica strain {j}, lag time = {str(tau_lag[1][j - 1])}')

        axs[0][s].set_ylim(-1.2, -0.4) # hard-coding ylim
        axs[0][s].set_ylabel('log(g + l)')
        axs[0][s].set_xlabel('Time')
        
        axs[1][s].plot(t_interval, locals()[f'M{s}'], label='methionine', color='m')
        axs[1][s].plot(t_interval, locals()[f'A{s}'], label='acetate', color='r')
        axs[1][s].plot(t_interval, locals()[f'L{s}'], label='lactose', color='y')

        axs[1][s].set_ylim(0, 20) #hard-coding ylim
        axs[1][s].set_ylabel('Amount of Nutrients')

        if s == 0: ## some hard code for titles and legends
            axs[0][s].set_title('Phase 1')
            axs[0][s].legend(loc='lower left')
            axs[1][s].legend(loc='lower left')
        else:
            axs[0][s].set_title('Phase 2')

    plt.show()

### run regular: observe as is
def regular(init_M, init_A, init_L, Ta, tau_lag=None, n=[1,1], frid=False):
    '''
    init_M: initial methionine concentration
    init_A: initial acetate concentration
    init_L: initial lactose concentration
    Ta: length of phase 1
    
    n:[nE, nS] where nE:int strains of E. coli, nS: "" for S. enterica
    tau_lag:[tau_lagE, tau_lagS] where tau_lagE:list of tau_lags for each strain of E. coli, tau_lagS: "" for S. enterica; or (default) lag time for each strain is a random value between 0 and 12 inclusive
    frid: when True, dMdt & dAdt & dLdt = 0 (for use with optimized())
    '''
    if tau_lag == None:
        tau_lag = [[random.uniform(0, 12) for i in range(n[0])], [random.uniform(0, 12) for j in range(n[1])]]
    else:
        n = [len(tau_lag[0]), len(tau_lag[1])]
        
    init_cond = [init_M, init_A, init_L] + [0, 1/(n[0] + n[1])]*(n[0] + n[1]) ## each strain for each species starts with growing population 0 and lagging population 1

    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, 15, 1000) #

    # run phases
    print(f'Running phase 1: antibiotic, with {init_cond}')
    sol1, para1 = run_phase(odes, init_cond, t_interval1, tau_lag, n)

    print(f'Running phase 2: no antibiotic, with {sol1[-1, :]}')
    sol2, para2 = run_phase(odes, sol1[-1, :], t_interval2, tau_lag, n, frid=frid, phase=2)

    #plot
    solplot([sol1, sol2], [para1, para2])
    
### run taus: Fridman analysis pt 1
def taus(init_M, init_A, init_L, Ta, n=[1,1], taulim=12, plot=True):
    '''
    init_M: initial methionine concentration
    init_A: initial acetate concentration
    init_L: initial lactose concentration
    Ta: length of phase 1

    n:[nE, nS] where nE:int strains of E. coli, nS: "" for S. enterica
    taulim: limit of tau_lags to try
    plot: when False, will not plot

    returns:list of max n*0, opt tau_lag
    '''
    tlim = 15 #

    # read in parameters
    nE = n[0]
    nS = n[1]

    taus = np.linspace(0, taulim, 100)[1:] ## lowered accuracy with 100 instead of 1000 for speed
    init_cond = [init_M, init_A, init_L] + [0, 1]*(n[0] + n[1]) ## each strain for each species starts with growing population 0 and lagging population 1

    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, tlim, 1000)

    ## run phases and append n*0 to list; use only one strain per species b/c all strains are identical in tau_lag
    n_0E = []
    n_0S = []
    
    for tau in taus:
        tau_lag = [[tau for i in range(1, nE + 1)], [tau for j in range(1, nS + 1)]]
        
        sol1, para1 = run_phase(odes, init_cond, t_interval1, tau_lag, n)
        sol2, para2 = run_phase(odes, sol1[-1, :], t_interval2, tau_lag, n, frid=True, phase=2)
        totalE = sol2[:, (1 + 2*1)][-1] + sol2[:, (2 + 2*1)][-1]
        totalS = sol2[:, ((1 + 2*nE) + 2*1)][-1] + sol2[:, ((2 + 2*nE) + 2*1)][-1]

        # append to n_0
        n_0E.append((math.e)**-tlim * totalE) ## calculate n*0 using the last t value available (== tlim)
        n_0S.append((math.e)**-tlim * totalS)
        

    # plot
    if plot:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        
        axs[0].plot(taus, n_0E, color = (0, 0.8, 0), label=f'E. coli strain 1')
        axs[1].plot(taus, n_0S, color = (0, 0, 0.8), label=f'S. enterica strain 1')

        for i in range(2):
            axs[i].set_ylabel('n*0')
            axs[i].set_xlabel('tau_lag')
            axs[i].legend()

        plt.show()
        
    return [[max(n_0E), taus[n_0E.index(max(n_0E))]], [max(n_0S), taus[n_0S.index(max(n_0S))]]]

### run 3d: Fridman analysis pt 2
def threed(init_M, init_A, init_L, n=[1,1], taulim=12, Talim=12, colors=None): ## colors is a temporary fix because its 10:27pm
    Tas = [Ta for Ta in range(1, Talim + 1)] # list of Tas to try, whole numbers only

    opt_tausE = []
    opt_tausS = []
    for Ta in Tas:
        opt_tausE.append(taus(init_M, init_A, init_L, Ta, n, taulim, plot=False)[0][1])
        opt_tausS.append(taus(init_M, init_A, init_L, Ta, n, taulim, plot=False)[1][1])

    # plot     
    plt.plot(Tas, opt_tausE, color = colors[0], label=f'E. init_M: {init_M} init_A: {init_A} init_L: {init_L}')
    plt.plot(Tas, opt_tausS, color = colors[1], label=f'S. init_M: {init_M} init_A: {init_A} init_L: {init_L}')

    plt.ylabel('Optimal tau_lag')
    plt.xlabel('Ta')
    
### user interface
inp = input("Run Sydney's test track (s) or run custom (c)?: ").strip()
if inp == 's':
    print('Running regular(10, 10, 10, 6, tau_lag=[[3, 6, 9], [6]], n=[3,1]). Initial methionine, initial acetate, initial lactose: 10; Ta: 6; tau_lag: 3,6,9 for E. coli and 6 for S. enterica; E. coli strains: 3; S. enterica strains: 1.')
    regular(10, 10, 10, 6, tau_lag=[[3, 6, 9], [6]], n=[3,1])
    print('Running taus(10, 10, 10, 6, n=[3,1]). Same parameters as above; Returns both max n*0 and optimal tau_lag for E. coli and S. enterica, respectively.')
    print(taus(10, 10, 10, 6, n=[3,1]))
    print('Running taus(10, 10, 10, 6, n=[5,1]).')
    print(taus(10, 10, 10, 6, n=[5,1]))
    print('Running taus(100, 100, 100, 6, n=[3,1]).')
    print(taus(100, 100, 100, 6, n=[3,1]))
    print('Running threed(10, 10, 10, n=[3,1], taulim=14) and threed(100, 100, 100, n=[3,1], taulim=12)')
    fig = plt.figure()
    ax = fig.add_subplot()
    threed(10, 10, 10, n=[3,1], taulim=14, colors=['lightgreen', 'lightblue'])
    threed(100, 100, 100, n=[3,1], taulim=12, colors=['darkgreen', 'royalblue'])
    ax.set_aspect('equal', adjustable='box')
    plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
    fig.tight_layout()
    plt.show()
    
elif inp == 'c':
    regular(10, 10, 10, 6, tau_lag=[[3, 6, 9], [6]], n=[3, 1])
    pass
