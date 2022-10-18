def odes_mono(x, t, alpha, lags, frid=False):
    nE = len(lags)

    # R
    M = x[0]
    L = x[1]

    # half-saturation constants
    K_M = 1
    K_L = 1

    # resource decay constants
    kM = 5e-9
    kL = 5e-9

    # resource constants
    cM = 0.1
    cL = 1.0

    # E growth constants
    rE = 1
    kE = 5e-9

    # E
    for i in range(nE):
        locals()[f'tau_lagE{i}'] = lags[i]

        locals()[f'El{i}'] = x[2 * i + 2]
        locals()[f'Eg{i}'] = x[2 * i + 3]

        # differential equations
        locals()[f'dEl{i}dt'] = -locals()[f'El{i}'] / locals()[f'tau_lagE{i}']
        locals()[f'dEg{i}dt'] = (1 - alpha) * rE * locals()[f'Eg{i}'] * (M / (M + K_M)) * (L / (L + K_L)) - kE * locals()[f'Eg{i}'] + locals()[f'El{i}'] / locals()[f'tau_lagE{i}']

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