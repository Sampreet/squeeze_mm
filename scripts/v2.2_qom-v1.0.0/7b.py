# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers import HLESolver, QCMSolver
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
            'min'   : 75.0,
            'max'   : 225.0,
            'dim'   : 301
        },
        'Y' : {
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
        'method'        : 'vode',
        'measure_codes' : ['entan_ln'],
        'indices'       : (0, 1),
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
        'x_label'           : '$n_{b}$',
        'x_tick_labels'     : ['$10^{' + str(i - 3) + '}$' for i in range(9)],
        'x_ticks'           : [10**(i - 3) for i in range(9)],
        'x_ticks_minor'     : sum([[10**(i - 3) * (j + 2) for i in range(8)] for j in range(7)], []),
        'x_scale'           : 'log',
        'y_colors'          : ['b', 'r', 'k'],
        'y_sizes'           : [2.0] * 2 + [1.0],
        'y_styles'          : ['-'] * 2 + [':'],
        'y_legend'          : [
            '$\\kappa = 0.1 \\omega_{m}$',
            '$\\kappa = 1.0 \\omega_{m}$'
        ],
        'v_label'           : '$E_{N_{\\mathrm{max}}}$',
        'v_tick_labels'     : ['{:0.1f}'.format(i * 0.1) for i in range(5)],
        'v_ticks'           : [i * 0.1 for i in range(5)],
        'v_ticks_minor'     : [i * 0.02 for i in range(21)],
        'show_legend'       : True,
        'legend_location'   : 'upper right',
        'height'            : 4.8,
        'width'             : 9.6,
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'annotations'       : [{
            'text'  : '(b)',
            'xy'    : (0.15, 0.84)
        }]
    }
}

# function to calculate the ratio and entanglement
def func_rat_entan_ln(system_params):
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

    # get modes, correlations and times
    Modes, Corrs, _ = HLESolver(
        system=system,
        params=params['solver']
    ).get_modes_corrs_times_in_range()
    # get entanglement
    eln = np.mean(QCMSolver(
        Modes=Modes,
        Corrs=Corrs,
        params=params['solver']
    ).get_measures(), axis=0)[0]

    return np.array([rat, eln])

if __name__ == '__main__':
    # low kappa
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/7b_kappa=0.1'
    params['system']['kappa_norm']          = 0.1
    looper  = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_entan_ln,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    xs      = looper.axes['Y']['val']
    elns_0  = np.max(looper.results['V'], axis=1).transpose()[1]

    # high kappa
    params['looper']['file_path_prefix']    = 'data/v2.2_qom-v1.0.0/7b_kappa=1.0'
    params['system']['kappa_norm']          = 1.0
    looper  = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func_rat_entan_ln,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    elns_1  = np.max(looper.results['V'], axis=1).transpose()[1]

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=xs,
        vs=[elns_0, elns_1]
    )
    plotter.show(
        hold=True
    )