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
        'show_progress' : True,
        'X' : {
            'var'   : 'beta_pm_sum',
            'min'   : 75,
            'max'   : 225,
            'dim'   : 301
        },
        'Y' : {
            'var'   : 'kappa_norm',
            'min'   : 1e-3,
            'max'   : 1e0,
            'dim'   : 151,
            'scale' : 'log'
        }
    },
    'solver': {
        'show_progress' : False,
        'cache'         : True,
        'method'        : 'vode',
        'indices'       : [(2, 2)],
        't_min'         : 0.0,
        't_max'         : 10000.0,
        't_dim'         : 100001,
        't_range_min'   : 99371,
        't_range_max'   : 100001
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
        'x_label'           : '$\\kappa / \\omega_{m}$',
        'x_ticks'           : [10**(i - 3) for i in range(4)],
        'x_ticks_minor'     : sum([[10**(i - 3) * (j + 2) for i in range(3)] for j in range(8)], []),
        'x_scale'           : 'log',
        'y_colors'          : ['b', 'r', 'k'],
        'y_legend'          : [
            '$n_{b} = 10$',
            '$n_{b} = 1000$'
        ],
        'y_sizes'           : [2.0] * 2 + [1.0],
        'y_styles'          : ['-'] * 2 + [':'],
        'y_name'            : '$n_{b}$',
        'v_label'           : '$\\frac{G_{1}}{G_{0}} |_{\\mathrm{opt}}$',
        'v_ticks'           : [0.6, 0.7, 0.8, 0.9, 1.0],
        'v_ticks_minor'     : [i * 0.02 + 0.6 for i in range(21)],
        'show_legend'       : True,
        'legend_location'   : 'upper right',
        'height'            : 4.8,
        'width'             : 9.6,
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'annotations'       : [{
            'text'  : '(b)',
            'xy'    : (0.18, 0.84)
        }]
    }
}

# function to calculate the ratio and variances
def func_rat_var(system_params):
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

    # get mechanical position variances
    var = np.min(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices_in_range()[:, 0])

    # update results
    return np.array([rat, var], dtype=np.float_)

if __name__ == '__main__':
    # low thermal phonons
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/5_n=10.0'
    params['system']['ns'][1]               = 10.0
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    xs          = looper.axes['Y']['val']
    _, _idxs_0  = np.argmin(looper.results['V'], axis=1).transpose()
    vs_0        = looper.results['V'].transpose()[0, _idxs_0, 0]

    # high thermal phonons
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/5_n=1000.0'
    params['system']['ns'][1]               = 1000.0
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    _, _idxs_1  = np.argmin(looper.results['V'], axis=1).transpose()
    vs_1        = looper.results['V'].transpose()[0, _idxs_1, 0]

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=xs,
        vs=[vs_0, vs_1]
    )
    plotter.show(
        hold=True
    )