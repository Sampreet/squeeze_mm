# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.MiddleMembrane import MM_01

# all parameters
params = {
    'looper': {
        'show_progress' : True,
        'X'             : {
            'var'   : 'beta_pm_sum',
            'min'   : 75.0,
            'max'   : 225.0,
            'dim'   : 301
        },
        'Y'             : {
            'var'   : 'ns',
            'idx'   : 1,
            'min'   : 1e-3,
            'max'   : 1e5,
            'dim'   : 81,
            'scale' : 'log'
        }
    },
    'solver': {
        'show_progress' : False,
        'cache'         : True,
        'ode_method'    : 'vode',
        'indices'       : [(2, 2)],
        't_min'         : 0.0,
        't_max'         : 1000.0,
        't_dim'         : 10001,
        't_index_min'   : 9371,
        't_index_max'   : 10001
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
        'colors'            : ['b', 'r', 'k'],
        'sizes'             : [2.0] * 2 + [1.0],
        'styles'            : ['-'] * 2 + [':'],
        'x_label'           : '$n_{b}$',
        'x_tick_labels'     : ['$10^{' + str(i - 3) + '}$' for i in range(9)],
        'x_ticks'           : [10**(i - 3) for i in range(9)],
        'x_ticks_minor'     : sum([[10**(i - 3) * (j + 2) for i in range(8)] for j in range(7)], []),
        'x_scale'           : 'log',
        'v_label'           : '$- 10 \\mathrm{log}_{10} \\left( \\langle \\tilde{Q}^{2} \\rangle \\right)$',
        'v_ticks'           : [i * 5 for i in range(5)],
        'v_ticks_minor'     : [i * 1 for i in range(21)],
        'show_legend'       : True,
        'legend_labels'     : [
            '$\\kappa = 0.1 \\omega_{m}$',
            '$\\kappa = 1.0 \\omega_{m}$'
        ],
        'legend_location'   : 'upper right',
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'width'             : 9.6,
        'height'            : 4.8,
        'annotations'   : [{
            'text'  : '(a)',
            'xy'    : (0.15, 0.84)
        }]
    }
}

# function to calculate the ratio and entanglement
def func_rat_entan_ln(system_params):
    # update parameters
    val = system_params['beta_pm_sum']
    system_params['betas'][1] = val / 2.0
    system_params['betas'][2] = val / 2.0

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

    # get mechanical position variance
    var = np.mean(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices(), axis=0)[0]

    return np.array([rat, var])

if __name__ == '__main__':
    # low kappa
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/7a_kappa=0.1'
    params['system']['kappa_norm'] = 0.1
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_entan_ln,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    xs  = looper.axes['Y']['val']
    vars_0 = np.min(looper.results['V'], axis=1).transpose()[1]

    # high kappa
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/7a_kappa=1.0'
    params['system']['kappa_norm'] = 1.0
    looper  = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_entan_ln,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    vars_1 = np.min(looper.results['V'], axis=1).transpose()[1]

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        vs=- 10 * np.log10([vars_0, vars_1]),
        xs=xs
    )
    plotter.show()