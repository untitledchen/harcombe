from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import random
import math

### define odes equation
def odes(x, t, alpha, nE, tau_lag, frid=False):
    '''
    x: list of initial conditions
    t: list of times
    alpha: additional death rate due to antibiotic, as a proportion of max growth rate r
    nE: number of strains of E. coli
    tau_lag: list of lag times of each strain
    frid: when True, dMdt and dLdt = 0 (for use with Fridman-style analysis)
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

        # solutions -- starting with Eg1 on x[2] and El1 on x[3], Eg and El occupy alternating indexes for each strain
        locals()[f'Eg{i}'] = x[2*i]
        locals()[f'El{i}'] = x[1 + 2*i]

        # differential equations
        locals()[f'dEg{i}dt'] = (1 - alpha)*locals()[f'rE{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - locals()[f'kE{i}']*locals()[f'Eg{i}'] + locals()[f'El{i}']/locals()[f'tau_lagE{i}']
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/locals()[f'tau_lagE{i}']

    # M
    sigma_cM = 0
    for i in range(1, nE + 1):
        sigma_cM += locals()[f'cM{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dMdt = (-sigma_cM - kM*M) * ([1, 0][frid])
    
    # L
    sigma_cL = 0
    for i in range(1, nE + 1):
        sigma_cL += locals()[f'cL{i}']*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L))
        
    dLdt = (-sigma_cL - kL*L) * ([1, 0][frid])

    to_return = [dMdt,dLdt]
    for i in range(1, nE + 1): ## could probably use list comprehension instead of this for loop
        to_return.append(locals()[f'dEg{i}dt'])
        to_return.append(locals()[f'dEl{i}dt'])

    return to_return
###

### run_phase 1 (antibiotic) or 2 (no antibiotic)
def run_phase(odes, init_cond, t_interval, nE, tau_lag, frid=False, phase=1):

    if phase == 1:
        alpha = 2
    elif phase == 2:
        alpha = 0
    
    sol = odeint(odes, init_cond, t_interval, args=(alpha, nE, tau_lag, frid)) ## when alpha = -2, effective growth rate = -r
    para = [t_interval, nE, tau_lag]
    return sol, para

### general matplotlib plot function
def solplot(sols, paras):
    '''
    sols:list of solutions for each phase
    paras:list of parameters for each phase [t_interval, nE, tau_lag]
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

        # convert absolute populations to frequencies of the total population, by creating totals and later freqs
        totals = locals()[f'Eg1{s}'] + locals()[f'El1{s}']
        for i in range(2, nE + 1):
            totals += locals()[f'Eg{i}{s}'] + locals()[f'El{i}{s}']
        
        g = 0.8/(1.5*nE)
        for i in range(1, nE + 1):
            freqs = (locals()[f'Eg{i}{s}']+locals()[f'El{i}{s}'])/totals
            axs[s].plot(t_interval, np.log10(freqs), color = (0, 0.8 - g*i, 0), label=f'E. coli strain {i}, lag time = {str(tau_lag[i - 1])}')

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
def regular(init_M, init_L, Ta, nE=1, tau_lag=None, frid=False):
    '''
    init_M: initial methionine concentration
    init_L: initial lactose concentration
    Ta: length of phase 1
    tau_lag:list of lag times for each strain; or (default) lag time for each strain is a random value between 0 and 12 inclusive
    frid: when True, dMdt and dLdt = 0 (for use with optimized())
    '''
    if tau_lag == None:
        tau_lag = [random.uniform(0, 12) for i in range(nE)]
    else:
        nE = len(tau_lag)
        
    init_cond = [init_M, init_L] + [0, 1]*nE ## each strain starts with growing population 0 and lagging population 1

    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, 15, 1000) #

    # run phases
    print(f'Running phase 1: antibiotic, with {init_cond}')
    sol1, para1 = run_phase(odes, init_cond, t_interval1, nE, tau_lag)

    print(f'Running phase 2: no antibiotic, with {sol1[-1, :]}')
    sol2, para2 = run_phase(odes, sol1[-1, :], t_interval2, nE, tau_lag, frid=frid, phase=2)

    #plot
    solplot([sol1, sol2], [para1, para2])
    
### run taus: Fridman analysis pt 1
def taus(init_M, init_L, Ta, nE=1, taulim=12, plot=True):
    '''
    init_M: initial methionine concentration
    init_L: initial lactose concentration
    Ta: length of phase 1
    
    taulim: limit of tau_lags to try
    plot: when False, will not plot

    returns:list of max n*0, opt tau_lag
    '''
    tlim = 15 #

    taus = np.linspace(0, taulim, 100)[1:] ## lowered accuracy with 100 instead of 1000 for speed
    init_cond = [init_M, init_L] + [0, 1]*nE ## each strain starts with growing population 0 and lagging population 1

    # time intervals
    t_interval1 = np.linspace(0, Ta, 1000)
    t_interval2 = np.linspace(0, tlim, 1000)

    # run phases and append n*0 to list(n_0)
    for i in range(1, nE + 1):
        locals()[f'n_0{i}'] = []
    
    for tau in taus:
        tau_lag = [tau for i in range(1, nE + 1)]
        
        sol1, para1 = run_phase(odes, init_cond, t_interval1, nE, tau_lag)
        sol2, para2 = run_phase(odes, sol1[-1, :], t_interval2, nE, tau_lag, frid=True, phase=2)
        for i in range(1, nE + 1):
            locals()[f'total{i}'] = sol2[:, (2*i)] + sol2[:, (1 + 2*i)]

        # append to n_0
        for i in range(1, nE + 1):
            locals()[f'n_0{i}'].append((math.e)**-tlim * (locals()[f'total{i}'][-1])) ## calculate n*0 using the last t value available (== tlim)

    # plot
    if plot:     
        g = 0.8/(2*nE)
        for i in range(1, nE + 1):
            plt.plot(taus, locals()[f'n_0{i}'], color = (0, 0.8 - g*i, 0), label=f'E. coli strain {i}, lag time = equal')

        plt.ylabel('n*0')
        plt.xlabel('tau_lag')
        plt.legend()

        plt.show()
        
    return [max(locals()[f'n_01']), taus[locals()[f'n_01'].index(max(locals()[f'n_01']))]]

### run 3d: Fridman analysis pt 2
def threed(init_M, init_L, nE=1, taulim=12, Talim=12):
    Tas = [Ta for Ta in range(1, Talim + 1)] # list of Tas to try, whole numbers only

    opt_taus = []
    for Ta in Tas:
        opt_taus.append(taus(init_M, init_L, Ta, nE, taulim, plot=False)[1])

    # plot
    plt.plot(Tas, opt_taus, label=f'init_M: {init_M} init_L: {init_L}')

    plt.ylabel('Optimal tau_lag')
    plt.xlabel('Ta')
    
### user interface
inp = input("Run Sydney's test track (s) or run custom (c)?: ").strip()
if inp == 's':
    print('Running regular(10, 10, 6, tau_lag=[3,6,9]). Initial methionine: 10; initial lactose: 10; Ta: 6; strains: 3; tau_lag: 3, 6, 9.')
    regular(10, 10, 6, tau_lag=[3,6,9])
    print('Running regular(10, 10, 6, nE=3). Initial methionine: 10; initial lactose: 10; Ta: 6; strains: 3; tau_lag: random.')
    regular(10, 10, 6, nE=3)
    print('Running regular(10, 10, 6, nE=3, frid=True). Same parameters but phase 2 is adapted for Fridman analysis.')
    regular(10, 10, 6, nE=3, frid=True)
    print('Running taus(10, 10, 6, nE=3). Returns max n*0 and optimal tau_lag, respectively.')
    print(taus(10, 10, 6, nE=3))
    print('Running taus(10, 10, 6, nE=5).')
    print(taus(10, 10, 6, nE=5))
    print('Running taus(100, 100, 6, nE=3).')
    print(taus(100, 100, 6, nE=3))
    print('Running threed(10, 10, nE=3, taulim=14) and threed(100, 100, nE=3, taulim=14)')
    fig = plt.figure()
    ax = fig.add_subplot()
    threed(10, 10, nE=3, taulim=14)
    threed(100, 100, nE=3, taulim=14)
    ax.set_aspect('equal', adjustable='box')
    plt.legend()
    plt.show()
    
elif inp == 'c':
    pass
