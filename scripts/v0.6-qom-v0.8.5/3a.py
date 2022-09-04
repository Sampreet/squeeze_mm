# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.utils.looper import wrap_looper
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'som-systems')))
# import system
from systems import MM_01

# all parameters
params = {
    'looper': {
        'show_progress_x': True,
        'X': {
            'var': 'beta_pm_sum',
            'min': 0,
            'max': 600,
            'dim': 601
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (2, 2),
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
        'x_label': '$G_{1} / G_{0}$',
        'x_ticks': [i * 0.1 + 0.5 for i in range(6)],
        'x_ticks_minor': [i * 0.025 + 0.5 for i in range(21)],
        'y_colors': ['b', 'r', 'k'],
        'y_sizes': [2.0] * 2 + [1.0],
        'y_styles': ['-'] * 2 + [':'],
        'y_name': '$n_{b}$',
        'v_label': '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks': [i * 6 for i in range(5)],
        'v_ticks_minor': [i * 3 for i in range(9)],
        'show_legend': True,
        'legend_location': 'upper center',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 28,
        'legend_font_size': 24,
        'tick_font_size': 24,
        'annotations': [{
            's': '(a)',
            'xy': (0.14, 0.84)
        }]
    }
}

# function to calculate the ratio and variance
def func_rat_var_min(system_params, val, logger, results):
    # update parameters
    system_params['betas'][1] = val / 3.5
    system_params['betas'][2] = val * 2.5 / 3.5
    system = MM_01(params=system_params)
    # extract parameters
    _, c = system.get_ivc()
    _len_D = 4 * system.num_modes**2
    _system_params = c[_len_D:] if len(c) > _len_D else c
    # get effective values
    G_0, _, G_p, _, _, _ = system.get_effective_values(params=_system_params)
    # get measure dynamics
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # update results
    results.append((val, [G_p / G_0, np.min(M)]))

# low thermal phonons
params['system']['ns'][1] = 10.0
looper = wrap_looper(SystemClass=None, params=params, func=func_rat_var_min, looper='XLooper', file_path_prefix='data/v0.6-qom-v0.8.5/3a_n=10')
rat, var_0 = np.transpose(looper.results['V']).tolist()

# high thermal phonons
params['system']['ns'][1] = 1000.0
looper = wrap_looper(SystemClass=None, params=params, func=func_rat_var_min, looper='XLooper', file_path_prefix='data/v0.6-qom-v0.8.5/3a_n=1000')
rat, var_1 = np.transpose(looper.results['V']).tolist()

# plotter
plotter = MPLPlotter(axes={
    'X': rat,
    'Y': [10, 1000]
}, params=params['plotter'])
plotter.update(xs=rat, vs=- 10 * np.log10([var_0, var_1, [0.5] * len(var_0)]))
plotter.show(True)