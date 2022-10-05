def odes(init_cond, t_interval, alpha, lags, frid=False):
    n = len(lags)
    
    M = init_cond[0]
    L = init_cond[1]

    # half-saturation constants
    K_M = 1
    K_L = 1

    # resource decay constants
    kM = 5e-9
    kL = 5e-9

    # growth constants
    r = 1 
    k = 5e-9

    # resource consumption/production constants
    cM = 0.1 
    cL = 1.0 

    # E
    for i in range(n):
        
        locals()[f'El{i}'] = init_cond[2*i]
        locals()[f'Eg{i}'] = init_cond[1 + 2*i]

        locals()[f'tau_lag{i}'] = lags[i]

        # differential equations
        locals()[f'dEg{i}dt'] =  (1 - alpha)*r*locals()[f'Eg{i}']*(M/(M+K_M))*(L/(L+K_L)) - k*locals()[f'Eg{i}'] + locals()[f'El{i}']/locals()[f'tau_lag{i}']
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}']/locals()[f'tau_lag{i}']

    # M
    sigma_strains = 0
    for i in range(n):
        sigma_strains += locals()[f'Eg{i}']
        
    dMdt = (-cM*sigma_strains*(M/(M+K_M))*(L/(L+K_L)) - kM*M) * ([1,0][frid])
    
    # L
    dLdt = (-cL*sigma_strains*(M/(M+K_M))*(L/(L+K_L)) - kL*L) * ([1,0][frid])

    to_return = [dMdt,dLdt]
    for i in range(n):
        to_return.append(locals()[f'dEl{i}dt'])
        to_return.append(locals()[f'dEg{i}dt'])

    return to_return
