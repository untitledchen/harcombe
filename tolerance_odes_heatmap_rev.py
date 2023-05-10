# M L A El1 Eg1 El2 Eg2 Sl1 Sg1
# x is a numpy array
# split x into resources as variables and populations g and l, as numpy arrays
import pdb

import numpy as np

def odes(x, t, alpha, lags, frid=False):
    #  resources vars
    M = x[0]
    L = x[1]
    A = x[2]

    # populations arrays
    nE = len(lags[0])

    ## assume np array!
    lagsE = np.array(lags[0])

    Egs = x[3:(nE + 3)]
    Els = x[(nE + 3):(2 * nE + 3)]

    alphaE = alpha[0]

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
    rE = 1
    kE = 5e-9
    # S
    rS = 0.5 ##
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

    dEls = -Els / lagsE
    dEgs = Egs * (1 - alphaE) * rE * (M/(M + K_M)) * (L/(L + K_L)) - Egs * kE + Els / lagsE

    # S. enterica
    if len(lags) > 1:
        nS = len(lags[1])
        lagsS = np.array(lags[1])
        alphaS = alpha[1]

        Sgs = x[(2 * nE + 3):(nS + 2 * nE + 3)]
        Sls = x[(nS + 2 * nE + 3):]

        dSls = -Sls / lagsS
        dSgs = Sgs * (1 - alphaS) * rS * (A / (A + K_A)) - Sgs * kS + Sls / lagsS
    else:
        Sgs = []
        Sls = []

        dSgs = []
        dSls = []

    # resource equations
    sigmaE = np.sum(Egs)
    sigmaS = np.sum(Sgs)

    # M
    dMdt = (-sigmaE * cM * (M / (M + K_M)) * (L / (L + K_L)) + sigmaS * pM * rS * (A/(A + K_A)) - kM * M) * (1, 0)[frid]

    # L
    dLdt = (-sigmaE * cL * (M / (M + K_M)) * (L / (L + K_L)) - kL * L) * (1, 0)[frid]

    # A
    dAdt = (sigmaE * pA * rE * (M / (M + K_M)) * (L / (L + K_L)) - sigmaS * cA * (A/(A + K_A)) - kA * A) * (1, 0)[frid]

    resources = [dMdt, dLdt, dAdt]
    to_return = np.concatenate((resources, dEgs, dEls, dSgs, dSls))

    return to_return

x = [1000, 1000, 0, 5, 0, 5, 0]
print(odes(x, 0, (3, 3), [[1], 1]))