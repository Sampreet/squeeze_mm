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
            'min'   : 75,
            'max'   : 225,
            'dim'   : 301
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
        'colors'            : ['b', 'r', 'b', 'r', 'k'],
        'sizes'             : [2.0] * 2 + [1.0] * 3,
        'styles'            : ['-'] * 2 + ['--'] * 2 + [':'],
        'x_label'           : '$G_{1} / G_{0}$',
        'x_ticks'           : [i * 0.1 + 0.5 for i in range(6)],
        'x_ticks_minor'     : [i * 0.01 + 0.5 for i in range(51)],
        'v_label'           : '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks'           : [i * 5 for i in range(5)],
        'v_ticks_minor'     : [i * 1 for i in range(21)],
        'show_legend'       : True,
        'legend_labels'     : [
            '$n_{b} = 10$',
            '$n_{b} = 1000$'
        ],
        'legend_location'   : 'upper center',
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'width'             : 9.6,
        'height'            : 4.8,
        'annotations'       : [{
            'text'  : '(a)',
            'xy'    : (0.15, 0.84)
        }]
    }
}

# function to calculate the ratio and variance
def func_rat_var(system_params):
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
    var = np.min(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices()[:, 0])

    # update results
    return np.array([rat, var], dtype=np.float_)

if __name__ == '__main__':
    # low thermal phonons with RWA
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/4a_rwa_n=10.0'
    params['system']['ns'][1] = 10.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    rats, vars_0_rwa = np.transpose(looper.results['V'])

    # high thermal phonons with RWA
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/4a_rwa_n=1000.0'
    params['system']['ns'][1] = 1000.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    _, vars_1_rwa = np.transpose(looper.results['V'])

    # low thermal phonons without RWA
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/4a_wrwa_n=10.0'
    params['system']['ns'][1] = 10.0
    params['system']['t_rwa'] = False
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    _, vars_0_wrwa = np.transpose(looper.results['V'])

    # high thermal phonons without RWA
    params['looper']['file_path_prefix'] = 'data/v2.2_qom-v1.0.1/4a_wrwa_n=1000.0'
    params['system']['ns'][1] = 1000.0
    params['system']['t_rwa'] = False
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_var,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    _, vars_1_wrwa = np.transpose(looper.results['V'])

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        vs=- 10 * np.log10([vars_0_rwa, vars_1_rwa, vars_0_wrwa, vars_1_wrwa, [0.5] * len(vars_0_rwa)]),
        xs=rats
    )
    plotter.show()