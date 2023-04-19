# M L A El1 Eg1 El2 Eg2 Sl1 Sg1
# x is a numpy array
# split x into resources as variables and populations g and l, as numpy arrays

def odes(x, t, alpha, lags, frid=False):
    #  resources vars
    M = x[0]
    L = x[1]
    A = x[2]

    # populations arrays
    n_species = len(lags)
    nE = len(lags[0])

    if n_species == 2:
        nS = len(lags[1])

        Sgs = x[(2 * nE + 3):(nS + 2 * nE + 3)]
        Sls = x[(nS + 2 * nE + 3):(2 * nS + 2 * nE + 3)]
    else:
        nS = 0

    Egs = x[3:(nE + 3)]
    Els = x[(nE + 3):(2*nE + 3)]

    ## CONSTANTS
    # half-saturation constants
    K_M = 1
    K_L = 1
    K_A = 1

    # resource decay constants
    kM = 5e-9
    kL = 5e-9
    kA = 5e-9

    # growth / decay constants
    # E
    alphaE = alpha[0]
    rE = 1
    kE = 5e-9
    # S
    rS = rs
    kS = 5e-9

    # resource consumption (c) / production (p) constants
    # E
    cM = 0.1
    cL = 1.0
    pA = 1.01
    # S
    cA = 1.0
    pM = 1.56

    ##



    for i in range(nE):
        locals()[f'tau_lagE{i}'] = lags[0][i]

        locals()[f'El{i}'] = x[2 * i + 3]
        locals()[f'Eg{i}'] = x[2 * i + 4]

        # differential equations
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}'] / locals()[f'tau_lagE{i}']
        locals()[f'dEg{i}dt'] = (1 - alphaE) * rE * locals()[f'Eg{i}'] * (M/(M + K_M)) * (L/(L + K_L)) - kE * locals()[f'Eg{i}'] + locals()[f'El{i}'] / locals()[f'tau_lagE{i}']

    # S. enterica
    if n_species == 2:
        alphaS = alpha[1]

    for j in range(nS):
        locals()[f'tau_lagS{j}'] = lags[1][j]

        locals()[f'Sl{j}'] = x[(3 + 2*nE) + 2*j]
        locals()[f'Sg{j}'] = x[(4 + 2*nE) + 2*j]

        locals()[f'dSl{j}dt'] = -locals()[f'Sl{j}'] / locals()[f'tau_lagS{j}']
        locals()[f'dSg{j}dt'] = (1 - alphaS)*rS*locals()[f'Sg{j}']*(A/(A + K_A)) - kS*locals()[f'Sg{j}'] + locals()[f'Sl{j}']/locals()[f'tau_lagS{j}']

    # resource equations
    sigma_strainsE = 0
    for i in range(nE):
        sigma_strainsE += locals()[f'Eg{i}']
    sigma_strainsS = 0
    for j in range(nS):
        sigma_strainsS += locals()[f'Sg{j}']

    # M
    dMdt = (-sigma_strainsE * cM * (M / (M + K_M)) * (L / (L + K_L)) + sigma_strainsS * pM * rS * (A/(A + K_A)) - kM * M) * (1, 0)[frid]

    # L
    dLdt = (-sigma_strainsE * cL * (M / (M + K_M)) * (L / (L + K_L)) - kL * L) * (1, 0)[frid]

    # A
    dAdt = (sigma_strainsE * pA * rE * (M / (M + K_M)) * (L / (L + K_L)) - sigma_strainsS * cA * (A/(A + K_A)) - kA * A) * (1, 0)[frid]

    to_return = [dMdt, dLdt, dAdt]
    for i in range(nE):
        to_return.append(locals()[f'dEl{i}dt'])
        to_return.append(locals()[f'dEg{i}dt'])
    for j in range(nS):
        to_return.append(locals()[f'dSl{j}dt'])
        to_return.append(locals()[f'dSg{j}dt'])

    return to_return