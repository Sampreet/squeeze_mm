# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper, run_loopers_in_parallel

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'som-systems')))
# import system
from systems import MM_01

# all parameters
params = {
    'looper': {
        'show_progress_y': True,
        'X': {
            'var': 'beta_pm_sum',
            'min': 0,
            'max': 600,
            'dim': 301
        },
        'Y': {
            'var': 'ns',
            'idx': 1,
            'min': -3,
            'max': 5,
            'dim': 201,
            'scale': 'log'
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'cache_dir': 'data/mm_01/0.0_1000.0_10001',
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
        'y_colors': ['b', 'r', 'k'],
        'y_sizes': [2.0] * 2 + [1.0],
        'y_styles': ['-'] * 2 + [':'],
        'y_name': '$\\kappa$',
        'y_unit': '$\\omega_{m}$',
        'v_label': '$E_{N_{\\mathrm{max}}}$',
        'v_tick_labels': ['{:0.1f}'.format(i * 0.1) for i in range(5)],
        'v_ticks': [i * 0.1 for i in range(5)],
        'v_ticks_minor': [i * 0.05 for i in range(9)],
        'show_legend': True,
        'legend_location': 'upper right',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 28,
        'legend_font_size': 24,
        'tick_font_size': 24,
        'annotations': [{
            's': '(b)',
            'xy': (0.13, 0.855)
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
    results.append((val, [G_p / G_0, np.max(M)]))

if __name__ == '__main__':
    # low kappa
    params['system']['kappa_norm'] = 0.1
    looper = wrap_looper(SystemClass=None, params=params, func=func_rat_var_min, looper='XYLooper', file_path_prefix='data/v1.4-qom-v0.9.0/7b_kappa=0.1')
    # looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_min, looper='XYLooper', file_path_prefix='data/v1.4-qom-v0.9.0/7b_kappa=0.1')
    ys = looper.axes['Y']['val']
    rat_0 = [np.transpose(v)[0][np.argmax(np.transpose(v)[1])] for v in looper.results['V']]
    eln_0 = [np.max(np.transpose(v)[1]) for v in looper.results['V']]

    # high kappa
    params['system']['kappa_norm'] = 1.0
    looper = wrap_looper(SystemClass=None, params=params, func=func_rat_var_min, looper='XYLooper', file_path_prefix='data/v1.4-qom-v0.9.0/7b_kappa=1.0')
    # looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_min, looper='XYLooper', file_path_prefix='data/v1.4-qom-v0.9.0/7b_kappa=1.0')
    rat_1 = [np.transpose(v)[0][np.argmax(np.transpose(v)[1])] for v in looper.results['V']]
    eln_1 = [np.max(np.transpose(v)[1]) for v in looper.results['V']]

    # plotter
    plotter = MPLPlotter(axes={
        'X': ys,
        'Y': [0.1, 1.0]
    }, params=params['plotter'])
    plotter.update(xs=ys, vs=[eln_0, eln_1])
    plotter.show(True)

    # # plotter
    # params['plotter']['v_label'] = '$\\frac{G_{1}}{G_{0}} |_{\\mathrm{opt}}$'
    # params['plotter']['v_tick_labels'] = [0.85, 0.9, 0.95, 1.0]
    # params['plotter']['v_ticks'] = [0.85, 0.9, 0.95, 1.0]
    # params['plotter']['v_ticks_minor'] = [i * 0.025 + 0.85 for i in range(7)]
    # plotter = MPLPlotter(axes={
    #     'X': ys,
    #     'Y': [0.1, 1.0]
    # }, params=params['plotter'])
    # plotter.update(xs=ys, vs=[rat_0, rat_1])
    # plotter.show(True)