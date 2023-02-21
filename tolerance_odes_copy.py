def odes_mono(x, t, alpha, lags, frid=False):
    nE = len(lags)

    # M = methionine, L = lactose
    M = x[0]
    L = x[1]

    # half-saturation constants
    K_M = 1
    K_L = 1

    # resource decay constants
    kM = 5e-9
    kL = 5e-9

    # resource consumption (c) / production (p) constants
    cM = 0.1
    cL = 1.0

    # E growth / decay constants
    alphaE = alpha[0]
    rE = 1
    kE = 5e-9

    # E
    for i in range(nE):
        locals()[f'tau_lagE{i}'] = lags[i]

        locals()[f'El{i}'] = x[2 * i + 2]
        locals()[f'Eg{i}'] = x[2 * i + 3]

        # differential equations
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}'] / locals()[f'tau_lagE{i}']
        locals()[f'dEg{i}dt'] = (1 - alphaE) * rE * locals()[f'Eg{i}'] * (M / (M + K_M)) * (L / (L + K_L)) - kE * locals()[f'Eg{i}'] + locals()[f'El{i}'] / locals()[f'tau_lagE{i}']

    # M
    sigma_strains = 0
    for i in range(nE):
        sigma_strains += locals()[f'Eg{i}']

    dMdt = (-sigma_strains * cM * (M / (M + K_M)) * (L / (L + K_L)) - kM * M) * ((1, 0)[frid])

    # L
    dLdt = (-sigma_strains * cL * (M / (M + K_M)) * (L / (L + K_L)) - kL * L) * ((1, 0)[frid])

    to_return = [dMdt, dLdt]
    for i in range(nE):
        to_return.append(locals()[f'dEl{i}dt'])
        to_return.append(locals()[f'dEg{i}dt'])

    return to_return

def odes_co(x, t, alpha, lags, frid=False):
    nE = len(lags[0])
    nS = len(lags[1])

    #  M = methionine, L = lactose, A = acetate
    M = x[0]
    L = x[1]
    A = x[2]

    # half-saturation constants
    K_M = 1
    K_L = 1
    K_A = 1

    # resource decay constants
    kM = 5e-9
    kL = 5e-9
    kA = 5e-9

    # E. coli
    # growth / decay constants
    alphaE = alpha[0]
    rE = 1
    kE = 5e-9

    # resource consumption (c) / production (p) constants
    cM = 0.1
    cL = 1.0
    pA = 1.01

    for i in range(nE):
        locals()[f'tau_lagE{i}'] = lags[0][i]

        locals()[f'El{i}'] = x[2 * i + 3]
        locals()[f'Eg{i}'] = x[2 * i + 4]

        # differential equations
        locals()[f'dEl{i}dt'] = round(-locals()[f'El{i}'] / locals()[f'tau_lagE{i}'], 5)
        locals()[f'dEg{i}dt'] = (1 - alphaE) * rE * locals()[f'Eg{i}'] * (M/(M + K_M)) * (L/(L + K_L)) - kE * locals()[f'Eg{i}'] + locals()[f'El{i}'] / locals()[f'tau_lagE{i}']

    # S. enterica
    # growth / decay constants
    alphaS = alpha[1]
    rS = 0.5
    kS = 5e-9

    # resource consumption (c) / production (p) constants
    cA = 1.0
    pM = 1.56

    for j in range(nS):
        locals()[f'tau_lagS{j}'] = lags[1][j]

        locals()[f'Sl{j}'] = x[(3 + 2*nE) + 2*j]
        locals()[f'Sg{j}'] = x[(4 + 2*nE) + 2*j]

        locals()[f'dSl{j}dt'] = round(-locals()[f'Sl{j}'] / locals()[f'tau_lagS{j}'], 5)
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