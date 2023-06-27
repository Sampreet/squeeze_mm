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
        'file_path_prefix'  : 'data/v2.2_qom-v1.0.0/6a',
        'X' : {
            'var'   : 'beta_pm_sum',
            'min'   : 75,
            'max'   : 225,
            'dim'   : 301
        }
    },
    'solver': {
        'show_progress' : False,
        'cache'         : False,
        'method'        : 'vode',
        'indices'       : [(2, 2)],
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
        'x_ticks_minor'     : [i * 0.01 + 0.5 for i in range(51)],
        'y_colors'          : ['b', 'r', 'k'],
        'y_sizes'           : [2.0] * 2 + [1.0],
        'y_styles'          : ['-'] * 2 + [':'],
        'y_legend'          : ['analytical', 'numerical'],
        'v_label'           : '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks'           : [i * 5 for i in range(5)],
        'v_ticks_minor'     : [i * 1 for i in range(21)],
        'show_legend'       : True,
        'legend_location'   : 'upper center',
        'height'            : 4.8,
        'width'             : 9.6,
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'annotations'       : [{
            'text'  : '(a)',
            'xy'    : (0.15, 0.84)
        }]
    }
}

# function to calculate the ratio and variance
def func_rat_vars(system_params):
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
    var = np.mean(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices_in_range()[:, 0])

    # get steady state variance
    var_ss = system.get_var_Q_ss_rwa(
        c=c
    )

    # get variance from FT
    var_ft = system.get_var_Q_ft_rwa(
        c=c
    )

    # update results
    return np.array([rat, var, var_ss, var_ft], dtype=np.float_)

# loop and plot
if __name__ == '__main__':
    # looper
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func_rat_vars,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    rats, vars, vars_ss, vars_ft = looper.results['V'].transpose()

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=rats,
        vs=- 10 * np.log10([vars_ss, vars_ft, [0.5] * len(vars_ss)])
    )
    plotter.add_scatter(
        xs=rats[::10],
        vs=- 10 * np.log10(vars[::10]),
        color='k',
        size=100,
        style='.'
    )
    plotter.show(
        hold=True
    )