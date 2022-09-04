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
            'var': 'alpha_pm',
            'min': 0.2,
            'max': 8.0,
            'dim': 40
        }
    },
    'solver': {
        'show_progress': False,
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
        'x_label': '$G_{1} / G_{0}$',
        'x_ticks': [0.5, 1],
        'y_colors': ['b', 'r'],
        'y_sizes': [2.0] * 2,
        'y_styles': ['-'] * 2,
        'y_name': '$n_{b}$',
        'v_label': '$- 10 \\mathrm{log}_{10} \\langle \\hat{Q}^{2} \\rangle$',
        'v_ticks': [i * 12 for i in range(3)],
        'show_legend': True,
        'legend_location': 'upper center',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 40,
        'legend_font_size': 32,
        'tick_font_size': 32,
        'annotations': [{
            's': '(a)',
            'xy': (0.18, 0.8)
        }]
    }
}

# function to calculate the minimum value of measure
def func_mdy_min(system_params, val, logger, results):
    # update parameters
    system_params['alphas'][1] = val
    system_params['alphas'][2] = val
    system = MM_01(params=system_params)
    alpha_0, alpha_m, alpha_p = system.params['alphas']
    beta_0, beta_m, beta_p = system.params['betas']
    val = (alpha_0 * (beta_m + beta_p) + 2 * alpha_p * beta_0) / (2 * alpha_0 * beta_0 + (alpha_m + alpha_p) * (beta_m + beta_p))
    # get measure dynamics
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # get minimum value
    m_min = np.min(M)
    # update results
    results.append((val, m_min))

# low thermal phonons
params['system']['ns'][1] = 10.0
looper = wrap_looper(SystemClass=None, params=params, func=func_mdy_min, looper='XLooper', file_path_prefix='data/v0.6-qom-v0.8.5/3a_alt_n=10')
xs = looper.results['X']
M_0 = looper.results['V']
print(np.max(- 10 * np.log10(M_0)))

# high thermal phonons
params['system']['ns'][1] = 1000.0
looper = wrap_looper(SystemClass=None, params=params, func=func_mdy_min, looper='XLooper', file_path_prefix='data/v0.6-qom-v0.8.5/3a_alt_n=1000')
M_1 = looper.results['V']
print(np.max(- 10 * np.log10(M_1)))

# plotter
plotter = MPLPlotter(axes={
    'X': xs,
    'Y': [10, 100],
}, params=params['plotter'])
plotter.update(xs=xs, vs=[(- 10 * np.log10(M_0)).tolist(), (- 10 * np.log10(M_1)).tolist()])
plotter.show(True)