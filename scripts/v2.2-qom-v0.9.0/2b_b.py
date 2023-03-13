# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'som-systems')))
# import system
from systems import MM_01

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': False,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (3, 3),
        'range_min': 19000,
        'range_max': 20001,
        't_min': 0.0,
        't_max': 200.0,
        't_dim': 20001
    },
    'system': {
        'alphas': [2.0, 0.2, 0.2],
        'betas': [100.0, 25.0, 25.0],
        'Delta_norm': 1.0,
        'g_norm': 1e-4,
        'gamma_norm': 1e-6,
        'kappa_norm': 0.1,
        'ns': [0.0, 10.0],
        'Omega_norms': [2.0, 2.0],
        't_rwa': True
    },
    'plotter': {
        'type': 'lines',
        'x_ticks': [190, 195, 200],
        'x_ticks_minor': [i * 1 + 190 for i in range(11)],
        'y_colors': ['b', 'r', 'k'],
        'y_sizes': [2.0] * 3,
        'y_styles': ['-', '-', '--'],
        'v_limits': [0.5, 1.1],
        'v_ticks': [0.6, 0.8, 1.0],
        'v_ticks_minor': [i * 0.05 + 0.5 for i in range(13)],
        'height': 3.2,
        'width': 4.0,
        'tick_font_size': 40
    }
}

# initialize logger
init_log()

# initialize system without RWA
params['system']['t_rwa'] = False
system = MM_01(params=params['system'])
# get measure dynamics
M_0, T = system.get_measure_dynamics(solver_params=params['solver'])

# initialize system with RWA
params['system']['t_rwa'] = True
system = MM_01(params=params['system'])
# get measure dynamics
M_1, T = system.get_measure_dynamics(solver_params=params['solver'])

# plotter
plotter = MPLPlotter(axes={
    'X': T,
    'Y': ['without RWA', 'with RWA', 'SQL']
}, params=params['plotter'])
plotter.update(xs=T, vs=[np.transpose(M_0)[0], np.transpose(M_1)[0], [0.5] * len(T), ])
plotter.show(True)