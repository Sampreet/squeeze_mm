# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers import HLESolver
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import run_loopers_in_parallel, wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.MiddleMembrane import MM_01

# all parameters
params = {
    'looper': {
        'show_progress'     : True,
        'X' : {
            'var'   : 'beta_pm_sum',
            'min'   : 75,
            'max'   : 225,
            'dim'   : 301
        }
    },
    'solver': {
        'show_progress' : False,
        'cache'         : True,
        'method'        : 'vode',
        'indices'       : [(2, 2), (3, 3)],
        't_min'         : 0.0,
        't_max'         : 1000.0,
        't_dim'         : 10001,
        't_range_min'   : 9371,
        't_range_max'   : 10001
    },
    'system': {
        'alphas'        : [2.0, 0.2, 0.2],
        'betas'         : [100.0, 25.0, 25.0],
        'Delta_norm'    : 1.0,
        'g_norm'        : 1e-4,
        'gamma_norm'    : 1e-6,
        'kappa_norm'    : 0.1,
        'ns'            : [0.0, 10.0],
        'Omega_norms'   : [2.0, 2.0],
        't_rwa'         : True
    },
    'plotter': {
        'type'              : 'lines',
        'x_label'           : '$G_{1} / G_{0}$',
        'x_ticks'           : [i * 0.1 + 0.5 for i in range(6)],
        'x_ticks_minor'     : [i * 0.02 + 0.5 for i in range(26)],
        'y_colors'          : ['b', 'r', 'k'],
        'y_sizes'           : [2.0] * 2 + [1.0],
        'y_styles'          : ['-', '-', ':'],
        'y_legend'          : [
            '$n_{b} = 10$',
            '$n_{b} = 1000$'
        ],
        'v_label'           : '$\\langle \\tilde{\\beta}^{\\dagger} \\tilde{\\beta} \\rangle$',
        'v_label_pad'       : -16,
        'v_scale'           : 'log',
        'v_tick_labels'     : ['$10^{' + str(i * 2 - 4) + '}$' for i in range(6)],
        'v_ticks'           : [10**(i * 2 - 4) for i in range(6)],
        'v_ticks_minor'     : [10.0**(i * 2 - 3) for i in range(11)],
        'show_legend'       : True,
        'legend_location'   : 'upper center',
        'height'            : 4.8,
        'width'             : 9.6,
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'annotations'       : [{
            'text'  : '(b)',
            'xy'    : (0.16, 0.84)
        }]
    }
}

# function to calculate the ratio and variance
def func_rat_n_beta(system_params):
    # update parameters
    val                         = system_params['beta_pm_sum']
    system_params['betas'][1]   = val / 2.0
    system_params['betas'][2]   = val / 2.0

    # initialize system
    system = MM_01(
        params=system_params
    )

    # get derived constants and controls
    _, _, c = system.get_ivc()
    
    # get squeezing ratio
    rat = system.get_params_ratio(
        c=c
    )

    # get mechanical position and momentum variances
    var_q, var_p = np.min(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices_in_range(), axis=0)

    # calculate hyperbolic angles
    r = np.arctanh(rat)
    chr = np.cosh(r)
    shr = np.sinh(r)

    # get phonon number in the Bogoluibov mode
    n_beta = (chr**2 + shr**2) * (var_q + var_p - 1) / 2.0 + shr**2 + chr * shr * (var_q - var_p)

    # update results
    return np.array([rat, n_beta], dtype=np.float_)

if __name__ == '__main__':
    # low thermal phonons
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/4b_n=10.0'
    params['system']['ns'][1]               = 10.0
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_n_beta,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    rats, n_betas_0 = np.transpose(looper.results['V'])

    # high thermal phonons
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/4b_n=1000.0'
    params['system']['ns'][1]               = 1000.0
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_n_beta,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    _, n_betas_1 = np.transpose(looper.results['V'])

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=rats,
        vs=[n_betas_0, n_betas_1, [1.0] * len(n_betas_0)]
    )
    plotter.show(
        hold=True
    )