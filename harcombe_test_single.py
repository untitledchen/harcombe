from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t, r1=1, K1=100, tau_lag1=6, q1=1, n1=1, r2=1, K2=100, tau_lag2=6, q2=1, n2=1):
    print(f'r1:{r1}/nK1:{K1}/ntau_lag1:{tau_lag1}/nq1:{q1}/nn1:{n1}/nr2:{r2}/nK2:{K2}/ntau_lag2:{tau_lag2}/nq2:{q2}/nn2:{n2}/n')
    ### E. coli
    R1 = x[0]

    for i in range(1, n1 + 1):

        # for now, set all r and K equal across strains
        locals()[f'r1{i}'] = r1
        locals()[f'K1{i}'] = K1
        locals()[f'tau_lag1{i}'] = tau_lag1
        locals()[f'q1{i}'] = q1

        locals()[f'N1{i}_g'] = x[2*i]
        locals()[f'N1{i}_l'] = x[1 + 2*i]

        interactions = 0
        for j in range(1, n1 + 1):
            if i == j:
                locals()[f'alpha1{i}{j}'] = 1
            else:
                locals()[f'alpha1{i}{j}'] = 1.5 ##

            interactions += locals()[f'alpha1{i}{j}']*locals()[f'N1{j}_g']
        
        locals()[f'dN1{i}_gdt'] = locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']) - interactions/locals()[f'K1{i}']) + locals()[f'N1{i}_l']/locals()[f'tau_lag1{i}']

        locals()[f'dN1{i}_ldt'] = - locals()[f'N1{i}_l']/locals()[f'tau_lag1{i}']
        print(locals())#

    ### S. Enterica
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
                locals()[f'alpha2{k}{m}'] = 1.5 ##

            interactions += locals()[f'alpha2{k}{m}']*locals()[f'N2{m}_g']
        
        locals()[f'dN2{k}_gdt'] = locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']) - interactions/locals()[f'K2{i}']) + locals()[f'N2{i}_l']/locals()[f'tau_lag2{i}']

        locals()[f'dN2{i}_ldt'] = - locals()[f'N2{i}_l']/locals()[f'tau_lag2{i}']

    dR1dt = 0
    for i in range(1, n1 + 1):
        dR1dt -= locals()[f'q1{i}']*locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']))

    for k in range(1, n2 + 1):
        dR1dt += locals()[f'q2{i}']*locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']))

    dR2dt = 0
    for i in range(1, n1 + 1):
        dR2dt += locals()[f'q1{i}']*locals()[f'r1{i}']*locals()[f'N1{i}_g']*(R1/(R1 + locals()[f'K1{i}']))

    for k in range(1, n2 + 1):
        dR2dt -= locals()[f'q2{i}']*locals()[f'r2{i}']*locals()[f'N2{i}_g']*(R2/(R2 + locals()[f'K2{i}']))

    to_return = [dR1dt, dR2dt]
    for i in range(1, n1 + 1):
        to_return.append(locals()[f'dN1{i}_gdt'])
        to_return.append(locals()[f'dN1{i}_ldt'])
    for k in range(1, n2 + 1):
        to_return.append(locals()[f'dN2{k}_gdt'])
        to_return.append(locals()[f'dN2{k}_ldt'])

    return to_return

    #dN1i_gdt = r1i*N1i_g*(R1/(R1 + K1i) - (N1i_g - alphaij*N1j_g)/K1i) + N1i_l/tau_lag1i
    #dR1dt = SUM-q1i*r1i*N1i_g*(R1/(R1 + K1i)) + SUM m*r2*N2_g*(R2/(R2 + K2))
    #dN2i_gdt = r2*N2_g*(R2/(R2 + K2) - (N2_g - alpha12*N1_g)/K1) + N1_l/tau_lag1

x0 = [10, 10, 0, 1, 0, 1]

print(odes(x=x0, t=0))
