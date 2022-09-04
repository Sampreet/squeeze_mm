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

# frequently used variables
_max = 3
_dim = 91

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': False,
        'method': 'zvode',
        'measure_type': 'wigner_s',
        'idx_e': [1],
        'xs': np.linspace(-_max, _max, _dim),
        'ys': np.linspace(-_max, _max, _dim),
        'range_min': 1980,
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
        'type': 'contourf',
        'x_label': '$Q$',
        'x_label_pad': -12,
        'x_tick_labels': [-_max, '', _max],
        'x_ticks': [-_max, 0, _max],
        'x_ticks_minor': [i - _max for i in range(7)],
        'y_label': '$P$',
        'y_label_pad': -12,
        'y_tick_labels': [-_max, '', _max],
        'y_ticks': [-_max, 0, _max],
        'y_ticks_minor': [i - _max for i in range(7)],
        'v_ticks': [0.0, 0.1, 0.2, 0.3],
        'v_ticks_minor': [i * 0.05 for i in range(8)],
        'height': 2.5,
        'width': 2.3,
        # 'show_cbar': True,
        # 'cbar_title': '$W$',
        # 'cbar_tick_pad': -2,
        # 'cbar_ticks': [0.0, 0.1, 0.2, 0.3],
        # 'width': 3.0,
    }
}

# initialize logger
init_log()

# initialize system without RWA
params['system']['t_rwa'] = False
system = MM_01(params=params['system'])
# get measure dynamics
M_0, T = system.get_measure_dynamics(solver_params=params['solver'])

for i in range(0, len(T), 5):
    params['plotter']['title'] = '$\\omega_{:} t = {:0.1f}$'.format('m', T[i])
    plotter = MPLPlotter(axes={
        'X': np.linspace(-_max, _max, _dim),
        'Y': np.linspace(-_max, _max, _dim)
    }, params=params['plotter'])
    plotter.update(vs=M_0[i][0])
    plotter.show(True)