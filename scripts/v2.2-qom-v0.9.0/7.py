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
            'min': 75.0,
            'max': 225.0,
            'dim': 301
        },
        'Y': {
            'var': 'ns',
            'idx': 1,
            'min': -3,
            'max': 5,
            'dim': 81,
            'scale': 'log'
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'cache_dir': 'H:/Workspace/data/mm_01/v2.2-qom-v0.9.0/7_0.0_1000.0_10001',
        'method': 'zvode',
        'range_min': 9900,
        'range_max': 10001,
        't_min': 0.0,
        't_max': 1000.0,
        't_dim': 10001
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
        'x_label': '$n_{b}$',
        'x_tick_labels': ['$10^{' + str(i - 3) + '}$' if i % 2 != 0 else '' for i in range(9)],
        'x_ticks': [10**(i - 3) for i in range(9)],
        'x_ticks_minor': sum([[10**(i - 3) * (j + 2) for i in range(8)] for j in range(7)], []),
        'x_scale': 'log',
        'y_colors': ['r', 'r', 'k', 'b', 'b'],
        'y_sizes': [2.0] * 2 + [1.0] + [2.0] * 2,
        'y_styles': ['-', '--'] + [':'] + ['-', '--'],
        'y_name': '$\\kappa$',
        'y_unit': '$\\omega_{m}$',
        'v_label_color': 'r',
        'v_label': '$- 10 \\mathrm{log}_{10} \\left( \\langle \\tilde{Q}^{2} \\rangle \\right)$',
        'v_ticks': [i * 5 for i in range(5)],
        'v_ticks_minor': [i * 1 for i in range(21)],
        'v_twin_label_color': 'b',
        'v_twin_label': '$E_{N_{\\mathrm{max}}}$',
        'v_twin_tick_labels': ['{:0.1f}'.format(i * 0.1) for i in range(5)],
        'v_twin_ticks': [i * 0.1 for i in range(5)],
        'v_twin_ticks_minor': [i * 0.02 for i in range(21)],
        'show_legend': True,
        'legend_location': 'upper right',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 32,
        'legend_font_size': 28,
        'tick_font_size': 28
    }
}

# function to calculate the ratio, variance and entanglement
def func_rat_var_entan(system_params, val, logger, results):
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
    params['solver']['measure_type'] = 'corr_ele'
    params['solver']['idx_e'] = (2, 2)
    m_var, _ = system.get_measure_dynamics(solver_params=params['solver'])
    params['solver']['measure_type'] = 'entan_ln'
    params['solver']['idx_e'] = (0, 1)
    m_entan, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # update results
    results.append((val, [G_p / G_0, np.min(m_var), np.max(m_entan)]))

if __name__ == '__main__':
    # low kappa
    params['system']['kappa_norm'] = 0.1
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_entan, looper='XYLooper', file_path_prefix='data/v2.2-qom-v0.9.0/7_kappa=0.1')
    xs = looper.axes['Y']['val']
    var_0 = [np.min(np.transpose(v)[1]) for v in looper.results['V']]
    entan_0 = [[np.max(np.transpose(v)[2])] for v in looper.results['V']]

    # high kappa
    params['system']['kappa_norm'] = 1.0
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_entan, looper='XYLooper', file_path_prefix='data/v2.2-qom-v0.9.0/7_kappa=1.0')
    var_1 = [np.min(np.transpose(v)[1]) for v in looper.results['V']]
    entan_1 = [[np.max(np.transpose(v)[2])] for v in looper.results['V']]

    # plotter
    plotter = MPLPlotter(axes={
        'X': xs,
        'Y': [0.1, 1.0]
    }, params=params['plotter'])
    plotter.update(xs=xs, vs=- 10 * np.log10([var_0, var_1, [0.5] * len(xs)]))
    plotter.update_twin_axis(xs=xs, vs=[entan_0, entan_1])
    plotter.show(True)