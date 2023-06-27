# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui import init_log
from qom.solvers import HLESolver
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.MiddleMembrane import MM_01

# frequently used variables
rat         = 0.93
beta_sum    = (400 * rat - 40.0) / (2.0 - 0.4 * rat)

# all parameters
params = {
    'solver': {
        'show_progress' : True,
        'cache'         : False,
        'method'        : 'vode',
        'indices'       : [(2, 2)],
        't_min'         : 0.0,
        't_max'         : 250.0,
        't_dim'         : 25001,
        't_range_min'   : 5000,
        't_range_max'   : 25001
    },
    'system': {
        'alphas'        : [2.0, 0.2, 0.2],
        'betas'         : [100.0, beta_sum / 2.0, beta_sum / 2.0],
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
        'x_label'           : '$\\omega_{m} t$',
        'x_ticks'           : [i * 50 + 50 for i in range(5)],
        'x_ticks_minor'     : [i * 5 + 50 for i in range(41)],
        'y_colors'          : ['k', 'k', 'k', 'b', 'r'],
        'y_sizes'           : [1.0] * 5,
        'y_styles'          : ['-.', '-', ':', '-', '-'],
        'v_label'           : '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks'           : [i * 5 for i in range(5)],
        'v_ticks_minor'     : [i * 1.25 for i in range(17)],
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

# initialize logger
init_log()

# initialize system without RWA
params['system']['t_rwa'] = False
system = MM_01(
    params=params['system']
)
# initialize solver
hle_solver  = HLESolver(
    system=system,
    params=params['solver']
)
# get times and mechanical position variances
T   = hle_solver.get_times_in_range()
M_0 = hle_solver.get_corr_indices_in_range().transpose()[0]

# initialize system with RWA
params['system']['t_rwa'] = True
system = MM_01(
    params=params['system']
)
# initialize solver
hle_solver  = HLESolver(
    system=system,
    params=params['solver']
)
# get mechanical position variances
M_1 = hle_solver.get_corr_indices_in_range().transpose()[0]

# get SQL
M_2 = [0.5] * len(T)

# frequently used variables
_, _, c = system.get_ivc()

# get analytical value
M_3 = [system.get_var_Q_ss_rwa(
    c=c
)] * len(T)

# get numerical value
M_4 = [system.get_var_Q_ft_rwa(
    c=c
)] * len(T)

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    xs=T,
    vs=- 10 * np.log10([M_0, M_1, M_2, M_3, M_4])
)
plotter.show(
    hold=True
)