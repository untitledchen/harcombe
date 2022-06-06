# locals() returns a dictionary of all variables in the local scope
# locals()['name'] = 0 creates a variable called 'name' with an assigned value of 0.

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t, r1=1, tau_lag1=6, r2=1, tau_lag2=6, K1=100, q1=1, n1=1, K2=100, q2=1, n2=1):
    '''
    x = [R1, R2,
        N11_g, N11_l... N1i_g, N1i_l,
        N21_g, N21_l... N2i_g, N2i_l]
    t: time interval
    r: growth rate
    tau_lag: lag time
    q: resource consumption constant
    n: number of strains
    '''

    ### 1: E. coli
    
    # R: resources
    R1 = x[0]

    # i: current strain
    for i in range(1, n1 + 1):

        ## for now, set all these constants equal across strains
        locals()[f'r1{i}'] = r1
        locals()[f'K1{i}'] = K1
        locals()[f'tau_lag1{i}'] = tau_lag1
        locals()[f'q1{i}'] = q1

        # N1i_g: growing population
        locals()[f'N1{i}_g'] = x[2*i] # occupies every other index in x starting from [2]
        
        # N1i_l: lagging population
        locals()[f'N1{i}_l'] = x[1 + 2*i] # occupies every other index in x starting from [3]
        
        # interactions with all other strains: + alpha1i2*N12 + alpha1i3*N13 + ... + alpha1ij*N1j
        ## would a matrix be better?

        interactions = 0
        for j in range(1, n1 + 1):
            if i == j:
                locals()[f'alpha1{i}{j}'] = 1 # when interacting with itself, alpha = 1 (derived from Lotka-Volterra Competition equations)
            else:
                locals()[f'alpha1{i}{j}'] = 1.5 ## interaction constant with strain j, arbitrary for now

            interactions += locals()[f'alpha1{i}{j}']*locals()[f'N1{j}_g'] ## wait why does this work? N1j_g should not exist at this point

        # dN1i_gdt
        locals()[f'dN1{i}_gdt'] = locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']) - interactions/locals()[f'K1{i}']) + locals()[f'N1{i}_l']/locals()[f'tau_lag1{i}']

        # dN1i_ldt
        locals()[f'dN1{i}_ldt'] = - locals()[f'N1{i}_l']/locals()[f'tau_lag1{i}']

    ### 2: S. Enterica
    ## identical to E. coli for now
    R2 = x[1]
    
    for k in range(1, n2 + 1):

        # for now, set all r and K equal across strains
        locals()[f'r2{i}'] = r2
        locals()[f'K2{i}'] = K2
        locals()[f'tau_lag2{i}'] = tau_lag2
        locals()[f'q2{i}'] = q2

        locals()[f'N2{k}_g'] = x[n1 + 2*k]
        locals()[f'N2{k}_l'] = x[n1 + 1 + 2*k]

        interactions = 0
        for m in range(1, n1 + 1):
            if k == m:
                locals()[f'alpha2{k}{m}'] = 1
            else:
                locals()[f'alpha2{k}{m}'] = 1.5

            interactions += locals()[f'alpha2{k}{m}']*locals()[f'N2{m}_g']
        
        locals()[f'dN2{k}_gdt'] = locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']) - interactions/locals()[f'K2{i}']) + locals()[f'N2{i}_l']/locals()[f'tau_lag2{i}']

        locals()[f'dN2{i}_ldt'] = - locals()[f'N2{i}_l']/locals()[f'tau_lag2{i}']

    # dR1dt: resources for E. coli
    dR1dt = 0
    # consumption from self and other strains
    for i in range(1, n1 + 1):
        dR1dt -= locals()[f'q1{i}']*locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']))
    # production from codependent species
    for k in range(1, n2 + 1):
        dR1dt += locals()[f'q2{i}']*locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']))

    # dR2dt: resources for S. enterica
    dR2dt = 0
    # production from codependent species
    for i in range(1, n1 + 1):
        dR2dt += locals()[f'q1{i}']*locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']))
    # consumption from self and other strains
    for k in range(1, n2 + 1):
        dR2dt -= locals()[f'q2{i}']*locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']))

    # create to_return in the form:
    #  [dR1dt, dR2dt,
    #  dN11_gdt, dN11_ldt... dN1i_gdt, dN1i_ldt,
    #  dN21_gdt, dN21_ldt... dN2i_gdt, dN2i_ldt]
    # similar to x from odes()
    to_return = [dR1dt, dR2dt]
    for i in range(1, n1 + 1):
        to_return.append(locals()[f'dN1{i}_gdt'])
        to_return.append(locals()[f'dN1{i}_ldt'])
    for k in range(1, n2 + 1):
        to_return.append(locals()[f'dN2{k}_gdt'])
        to_return.append(locals()[f'dN2{k}_ldt'])

    return to_return

#print(odes(x=x0, t=0))

### n*0 (approximation) over tau_lag
n_0 = [] # list of all calculated n*0s for each tau_lag tested
##taus = np.linspace(0, 12, 1000)[1:] # list of tau_lags to be tested, ranging from just above 0 to 12.00
    # [1:] to avoid a divide by zero error
taus = [6]

for tau in taus:
    # phase 1, antibiotics
    x0 = [10, 10, 0, 1, 0, 1]
    t_interval = np.linspace(0, 6, 1000)
    x = odeint(odes, x0, t_interval, args=(-1, tau, -1, tau))

    R1 = x[:, 0]
    R2 = x[:, 1]
    G1 = x[:, 2]
    L1 = x[:, 3]
    G2 = x[:, 4]
    L2 = x[:, 5]
    '''
    plt.semilogy(t_interval, R1, label='R1')
    plt.semilogy(t_interval, R2, label='R2')
    
    plt.semilogy(t_interval, G1, label='G1')
    plt.semilogy(t_interval, L1, label='L1')
    plt.semilogy(t_interval, G2, label='G2')
    plt.semilogy(t_interval, L2, label='L2')
    '''
    
    plt.semilogy(t_interval, G1+L1, label='1')
    plt.semilogy(t_interval, G2+L2, label='2')

    plt.title('Antibiotic Phase, a = -1')
    plt.ylabel('log(g+l)')
    plt.xlabel('t')
    plt.legend()

    plt.show()
'''
    # phase 2, no antibiotics
    x2 = [R1[-1], R2[-1], G1[-1], L1[-1], G2[-1], L2[-1]]
    t_interval2 = np.linspace(0, 15, 1000)
    y = odeint(odes, x2, t_interval2, args=(1, tau, 1, tau))

    R1 = x[:, 0]
    R2 = x[:, 1]
    G1 = x[:, 2]
    L1 = x[:, 3]
    G2 = x[:, 4]
    L2 = x[:, 5]

    plt.semilogy(t_interval, G1, label='G1')
    plt.semilogy(t_interval, L1, label='L1')
    plt.semilogy(t_interval, G2, label='G2')
    plt.semilogy(t_interval, L2, label='L2')

    plt.title('Non-antibiotic Phase, a = 1')
    plt.ylabel('log(g+l)')
    plt.xlabel('t')
    plt.ylim(sum(x2)*1e-2) # set the y-axis range to a reasonable distance below 0, so we can see where the asymptotic line would approximately intersect the y-axis
    plt.legend()

    plt.show()
'''
    ##total = G1 + L1

    # append calculated n*0 to n_0
    ##n_0.append((math.e)**-tlim2 * (total[-1])) # calculate n*0 using the last t value available (== tlim2)
'''
# plot: tau_lag vs. calculated n*0
plt.plot(taus, n_0)

plt.title('Effective Population (n*0) as a Result of Lag Time (tau_lag)')
plt.ylabel('n*0')
plt.xlabel('tau_lag')

plt.show()
'''
    # base versions of our equations, kept down here for my reference
    #dN1i_gdt = r1i*N1i_g*(R1/(R1 + K1i) - (N1i_g - alphaij*N1j_g)/K1i) + N1i_l/tau_lag1i
    #dR1dt = SUM-q1i*r1i*N1i_g*(R1/(R1 + K1i)) + SUM m*r2*N2_g*(R2/(R2 + K2))
    #dN2i_gdt = r2*N2_g*(R2/(R2 + K2) - (N2_g - alpha12*N1_g)/K1) + N1_l/tau_lag1
