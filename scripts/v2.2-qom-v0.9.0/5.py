# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import run_loopers_in_parallel

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
            'min': 75,
            'max': 225,
            'dim': 301
        },
        'Y': {
            'var': 'kappa_norm',
            'min': -3,
            'max': 0,
            'dim': 91,
            'scale': 'log'
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'cache_dir': 'H:/Workspace/data/mm_01/v2.2-qom-v0.9.0/5_0.0_10000.0_100001',
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (2, 2),
        'range_min': 99900,
        'range_max': 100001,
        't_min': 0.0,
        't_max': 10000.0,
        't_dim': 100001
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
        'x_label': '$\\kappa / \\omega_{m}$',
        'x_ticks': [10**(i - 2) for i in range(3)],
        'x_ticks_minor': sum([[10**(i - 2) * (j + 2) for i in range(2)] for j in range(8)], []),
        'x_scale': 'log',
        'y_colors': ['r', 'r', 'k', 'b', 'b'],
        'y_sizes': [2.0] * 2 + [1.0] + [1.0] * 2,
        'y_styles': ['-', '--'] + [':'] + ['-', '--'],
        'y_name': '$n_{b}$',
        'y_legend': ['$n_{b} = 10$', '$n_{b} = 1000$'],
        'v_label_color': 'r',
        'v_label': '$- 10 \\mathrm{log}_{10} \\left( \\langle \\tilde{Q}^{2} \\rangle \\right)$',
        'v_limits': [3, 22],
        'v_ticks': [i * 5 + 5 for i in range(4)],
        'v_ticks_minor': [i * 1 + 3 for i in range(20)],
        'v_twin_label_color': 'b',
        'v_twin_label': '$\\frac{G_{1}}{G_{0}} |_{\\mathrm{opt}}$',
        'v_twin_limits': [0.66, 1.04],
        'v_twin_ticks': [0.7, 0.8, 0.9, 1.0],
        'v_twin_ticks_minor': [i * 0.02 + 0.66 for i in range(20)],
        'show_legend': True,
        'legend_location': 'upper right',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 32,
        'legend_font_size': 28,
        'tick_font_size': 28
    }
}

# function to calculate the ratio and variance
def func_rat_var(system_params, val, logger, results):
    # update parameters
    system_params['betas'][1] = val / 2.0
    system_params['betas'][2] = val / 2.0
    system = MM_01(params=system_params)
    # extract parameters
    _, c = system.get_ivc()
    _len_D = 4 * system.num_modes**2
    _system_params = c[_len_D:] if len(c) > _len_D else c
    # get effective values
    G_0, _, G_p, _, _, _ = system.get_effective_values(params=_system_params)
    # get measure ft
    m_ft = system.get_var_Q_ft_rwa(params=_system_params)
    if m_ft != m_ft:
        m_ft = 100
    # get measure dynamics
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    m_min = np.min([100 if val != val else val for val in np.transpose(M)[0]])
    # get measure ss
    m_ss = system.get_var_Q_ss_rwa(params=_system_params)
    if m_ss != m_ss:
        m_ss = 100
    # update results
    results.append((val, [G_p / G_0, m_ft, m_min, m_ss]))

if __name__ == '__main__':
    # low thermal phonons with rwa
    params['system']['ns'][1] = 10.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var, looper='XYLooper', file_path_prefix='data/v2.2-qom-v0.9.0/5_rwa_n=10')
    xs = looper.axes['Y']['val']
    vs_0_var = list()
    vs_0_opt = list()
    for V in looper.results['V']:
        rat, m_ft, m_min, m_ss = np.transpose(V).tolist()
        vs_0_var.append(np.min(m_min))
        vs_0_opt.append(rat[np.argmin(m_min)])

    # high thermal phonons with rwa
    params['system']['ns'][1] = 1000.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var, looper='XYLooper', file_path_prefix='data/v2.2-qom-v0.9.0/5_rwa_n=1000')
    xs = looper.axes['Y']['val']
    vs_1_var = list()
    vs_1_opt = list()
    for V in looper.results['V']:
        _, m_ft, m_min, m_ss = np.transpose(V).tolist()
        vs_1_var.append(np.min(m_min))
        vs_1_opt.append(rat[np.argmin(m_min)])

    # plotter
    plotter = MPLPlotter(axes={
        'X': xs
    }, params=params['plotter'])
    plotter.update(xs=xs, vs=- 10 * np.log10([vs_0_var, vs_1_var, [0.5] * len(vs_0_var)]))
    plotter.update_twin_axis(xs=xs, vs=[vs_0_opt, vs_1_opt])
    plotter.show(True)