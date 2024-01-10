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
        'indices'       : [(3, 3)],
        't_min'         : 0.0,
        't_max'         : 200.0,
        't_dim'         : 2001,
        't_index_min'   : 1900,
        't_index_max'   : 2001
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
        'type'          : 'lines',
        'colors'        : ['b', 'r', 'k'],
        'sizes'         : [2.0] * 3,
        'styles'        : ['-', '-', '--'],
        'x_ticks'       : [190, 195, 200],
        'x_ticks_minor' : [i * 1 + 190 for i in range(11)],
        'v_limits'      : [0.5, 1.1],
        'v_ticks'       : [0.6, 0.8, 1.0],
        'v_ticks_minor' : [i * 0.05 + 0.5 for i in range(13)],
        'width'         : 4.0,
        'height'        : 3.2,
        'tick_font_size': 40
    }
}

# initialize logger
init_log()

# initialize system without RWA
params['system']['t_rwa'] = False
system = MM_01(
    params=params['system']
)

# get mechanical momentum variances
hle_solver = HLESolver(
    system=system,
    params=params['solver']
)
T = hle_solver.get_times()
M_0 = hle_solver.get_corr_indices().transpose()[0]

# initialize system with RWA
params['system']['t_rwa'] = True
system = MM_01(
    params=params['system']
)

# get mechanical momentum variances
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