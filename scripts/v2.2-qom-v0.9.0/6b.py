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

rat = 0.93
beta_sum = (400 * rat - 40.0) / (2.0 - 0.4 * rat)
print(beta_sum)
# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': False,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (2, 2),
        'range_min': 0,
        'range_max': 100001,
        't_min': 0.0,
        't_max': 1000.0,
        't_dim': 100001
    },
    'system': {
        'alphas': [2.0, 0.2, 0.2],
        'betas': [100.0, beta_sum / 2.0, beta_sum / 2.0],
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
        'x_ticks': [i * 50 + 50 for i in range(5)],
        'x_ticks_minor': [i * 5 + 50 for i in range(41)],
        'y_colors': ['k', 'k', 'k', 'b', 'r'],
        'y_sizes': [1.0] * 5,
        'y_styles': ['-.', '-', ':', '-', '-'],
        'v_label': '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks': [i * 5 for i in range(5)],
        'v_ticks_minor': [i * 1.25 for i in range(17)],
        # 'show_legend': True,
        # 'legend_location': 'upper right',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 32,
        'legend_font_size': 28,
        'tick_font_size': 28,
        'annotations': [{
            's': '(b)',
            'xy': (0.15, 0.84)
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
    'Y': ['without RWA', 'with RWA', 'SQL']
}, params=params['plotter'])
plotter.update(xs=T, vs=- 10 * np.log10([np.transpose(M_0)[0], np.transpose(M_1)[0], [0.5] * len(T), [0.031429448904571625] * len(T), [0.026700674850321043] * len(T)]))
plotter.show(True)