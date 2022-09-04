# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'som-systems')))
# import system
from systems import MM_01

# all parameters
params = {
    'looper': {
        'show_progress_x': True,
        'X': {
            'var': 'ns',
            'idx': 1,
            'min': -3,
            'max': 5,
            'dim': 401,
            'scale': 'log'
        },
        'Y': {
            'var': 'kappa_norm',
            'val': [0.1, 1.0]
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'method': 'zvode',
        'measure_type': 'entan_ln',
        'idx_e': (0, 1),
        'range_min': 9900,
        'range_max': 10001,
        't_min': 0.0,
        't_max': 1000.0,
        't_dim': 10001
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
        'x_label': '$n_{b}$',
        'x_tick_labels': ['$10^{' + str(i - 3) + '}$' for i in range(9)],
        'x_ticks': [10**(i - 3) for i in range(9)],
        'x_ticks_minor': sum([[10**(i - 3) * (j + 2) for i in range(8)] for j in range(7)], []),
        'x_scale': 'log',
        'y_colors': ['b', 'r'],
        'y_sizes': [2.0] * 2,
        'y_styles': ['-'] * 2,
        'y_name': '$\\kappa$',
        'y_unit': '$\\omega_{m}$',
        'v_label': '$E_{N} \\times 10^{3}$',
        'v_tick_labels': [i * 3 for i in range(5)],
        'v_ticks': [i * 0.003 for i in range(5)],
        'v_ticks_minor': [i * 0.0015 for i in range(9)],
        'show_legend': True,
        'legend_location': 'upper center',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 28,
        'legend_font_size': 24,
        'tick_font_size': 24,
        'annotations': [{
            's': '(b)',
            'xy': (0.14, 0.84)
        }]
    }
}

# looper
looper = wrap_looper(SystemClass=MM_01, params=params, func='mav', looper='XYLooper', file_path_prefix='data/v0.6-qom-v0.8.5/6b')
xs = looper.axes['X']['val']
vs = looper.results['V']

# plotter
plotter = MPLPlotter(axes={
    'X': xs,
    'Y': [0.1, 1.0]
}, params=params['plotter'])
plotter.update(xs=xs, vs=[vs[0], vs[1]])
plotter.show(True)