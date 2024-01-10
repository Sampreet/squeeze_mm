# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.solvers.measure import get_Wigner_distributions_single_mode
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.MiddleMembrane import MM_01

# frequently used variables
_max = 3
_dim = 601

# all parameters
params = {
    'solver': {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        'indices'       : [1],
        'wigner_xs'     : np.linspace(-_max, _max, _dim),
        'wigner_ys'     : np.linspace(-_max, _max, _dim),
        't_min'         : 0.0,
        't_max'         : 200.0,
        't_dim'         : 2001,
        't_index_min'   : 1980,
        't_index_max'   : 2001,
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
        'type'          : 'contourf',
        'x_label'       : '$Q$',
        'x_label_pad'   : -12,
        'x_tick_labels' : [-_max, '', _max],
        'x_ticks'       : [-_max, 0, _max],
        'x_ticks_minor' : [i - _max for i in range(7)],
        'y_label'       : '$P$',
        'y_label_pad'   : -12,
        'y_tick_labels' : [-_max, '', _max],
        'y_ticks'       : [-_max, 0, _max],
        'y_ticks_minor' : [i - _max for i in range(7)],
        'v_ticks'       : [0.0, 0.1, 0.2, 0.3],
        'v_ticks_minor' : [i * 0.05 for i in range(8)],
        'show_cbar'     : True,
        'cbar_title'    : '$W$',
        'cbar_tick_pad' : -2,
        'cbar_ticks'    : [0.0, 0.1, 0.2, 0.3],
        'width'         : 2.75,
        'height'        : 2.5
    }
}

# initialize logger
init_log()

# initialize system with RWA
params['system']['t_rwa'] = True
system = MM_01(
    params=params['system']
)

# get times and correlations
hle_solver = HLESolver(
    system=system,
    params=params['solver']
)
T = hle_solver.get_times()
_, Corrs = hle_solver.get_modes_corrs()
# get Wigner distributions
Wigners = get_Wigner_distributions_single_mode(
    Corrs=Corrs,
    params=params['solver']
)

# plotter
for i in range(0, len(T), 5):
    params['plotter']['title'] = '$\\omega_{:} t = {:0.1f}$'.format('m', T[i])
    plotter = MPLPlotter(
        axes={
            'X': np.linspace(-_max, _max, _dim),
            'Y': np.linspace(-_max, _max, _dim)
        },
        params=params['plotter']
    )
    plotter.update(
        vs=Wigners[i, 0]
    )
    plotter.show()