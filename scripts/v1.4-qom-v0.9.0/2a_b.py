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
        'idx_e': (2, 2),
        'range_min': 0,
        'range_max': 2001,
        't_min': 0.0,
        't_max': 200.0,
        't_dim': 2001
    },
    'system': {
        'alphas': [2.0, 0.8, 0.8],
        'betas': [100.0, 25.0, 62.5],
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
        'x_label': '$\\omega_{m} t$',
        'x_tick_labels': [0, 130, 140, 150, 160],
        'x_ticks': [i * 10 + 120 for i in range(5)],
        'x_ticks_minor': [i * 2 + 120 for i in range(21)],
        'y_colors': ['k', 'b', 'r'],
        'y_sizes': [1.0] * 3,
        'y_styles': ['--', '-', '-'],
        'v_label': '$\\langle \\hat{Q}^{2} \\rangle$',
        'v_ticks': [i * 2 for i in range(3)],
        'v_ticks_minor': [i * 0.5 for i in range(9)],
        'show_legend': True,
        'legend_location': 'upper right',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 40,
        'legend_font_size': 32,
        'tick_font_size': 32,
        'annotations': [{
            's': '(a)',
            'xy': (0.16, 0.8)
        }]
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
    'Y': ['SQL', 'without RWA', 'with RWA']
}, params=params['plotter'])
plotter.update(xs=T, vs=[[0.5] * len(T), np.transpose(M_0)[0], np.transpose(M_1)[0]])
plotter.show(True)