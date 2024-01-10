# dependencies
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.MiddleMembrane import MM_01

# all parameters
params = {
    'solver': {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        'indices'       : [(2, 2)],
        't_min'         : 0.0,
        't_max'         : 200.0,
        't_dim'         : 2001
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
        'sizes'             : [1.0] * 3,
        'styles'            : ['-', '-', '--'],
        'x_label'           : '$\\omega_{m} t$',
        'x_ticks'           : [i * 40 for i in range(6)],
        'x_ticks_minor'     : [i * 5 for i in range(41)],
        'v_label'           : '$\\langle \\hat{Q}^{2} \\rangle$',
        'v_ticks'           : [i * 1 for i in range(5)],
        'v_ticks_minor'     : [i * 0.25 for i in range(17)],
        'show_legend'       : True,
        'legend_labels'     : ['without RWA', 'with RWA', 'SQL'],
        'legend_location'   : 'upper right',
        'label_font_size'   : 32,
        'legend_font_size'  : 28,
        'tick_font_size'    : 28,
        'width'             : 9.6,
        'height'            : 4.8,
        'annotations'       : [{
            'text'  : '(a)',
            'xy'    : (0.14, 0.82)
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
hle_solver = HLESolver(
    system=system,
    params=params['solver']
)
# get times and variances
T = hle_solver.get_times()
M_0 = hle_solver.get_corr_indices().transpose()[0]

# initialize system with RWA
params['system']['t_rwa'] = True
system = MM_01(
    params=params['system']
)

# get variances
M_1 = HLESolver(
    system=system,
    params=params['solver']
).get_corr_indices().transpose()[0]

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    vs=[M_0, M_1, [0.5] * len(T)],
    xs=T
)
plotter.show()