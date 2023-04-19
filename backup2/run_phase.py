import numpy as np
from scipy.integrate import odeint
from backup2.tolerance_odes_copy import odes


def run_phase(alpha, init_cond, lags, t, phase, inc=1000, frid=False): #
    alpha_this = tuple([[a,0][phase-1] for a in alpha])
    t_interval = np.linspace(0, t, inc)
    sol = odeint(odes, init_cond, t_interval, args=(alpha_this, lags, frid))

    # plt.plot(t_interval, sol[:, 0], label='met', color='r')
    # plt.plot(t_interval, sol[:, 1], label='lac', color='b')
    # plt.plot(t_interval, sol[:, 2], label='ace', color='y')
    # plt.legend(loc='center right')
    # plt.title(phase)
    # plt.show()

    return sol