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
            'min': 75,
            'max': 225,
            'dim': 601
        }
    },
    'solver': {
        'show_progress': False,
        'cache': True,
        'cache_dir': 'H:/Workspace/data/mm_01/v2.2-qom-v0.9.0/4_0.0_1000.0_10001',
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
        'x_label': '$G_{1} / G_{0}$',
        'x_ticks': [i * 0.1 + 0.5 for i in range(6)],
        'x_ticks_minor': [i * 0.01 + 0.5 for i in range(51)],
        'y_colors': ['b', 'r', 'k'],
        'y_sizes': [2.0] * 2 + [1.0],
        'y_styles': ['-'] * 2 + [':'],
        'y_legend': ['analytical', 'numerical'],
        'v_label': '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_ticks': [i * 5 for i in range(5)],
        'v_ticks_minor': [i * 1 for i in range(21)],
        'show_legend': True,
        'legend_location': 'upper center',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 32,
        'legend_font_size': 28,
        'tick_font_size': 28,
        'annotations': [{
            's': '(a)',
            'xy': (0.15, 0.84)
        }]
    }
}

# function to calculate the ratio and variance
def func_rat_ss_ft(system_params, val, logger, results):
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
    # get steady-state variance
    var_ss = system.get_var_Q_ss_rwa(params=_system_params)
    var_ft = system.get_var_Q_ft_rwa(params=_system_params)
    # update results
    results.append((val, [G_p / G_0, var_ss, var_ft]))

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
    rat = G_p / G_0
    # get measure dynamics
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # variance
    var = np.min(np.transpose(M)[0])
    # update results
    results.append((val, [rat, var]))

if __name__ == '__main__':
    # low thermal phonons analytical and numerical
    looper = wrap_looper(SystemClass=None, params=params, func=func_rat_ss_ft, looper='XLooper')
    rat, var_ss, var_ft = np.transpose(looper.results['V'])

    # low thermal phonons with RWA
    looper = wrap_looper(SystemClass=None, params=params, func=func_rat_var, looper='XLooper', file_path_prefix='data/v2.2-qom-v0.9.0/4_rwa_n=10')
    rat, var, _ = np.transpose(looper.results['V']).tolist()

    # plotter
    plotter = MPLPlotter(axes={
        'X': rat,
        'Y': [0, 1]
    }, params=params['plotter'])
    plotter.update(xs=rat, vs=- 10 * np.log10([var_ss, var_ft, [0.5] * len(var_ss)]))
    ax = plotter.get_current_axis()
    ax.scatter(rat[::10], - 10 * np.log10(var)[::10], c='k', s=50)
    plotter.show(True)
    # print(var_ss[258], var_ft[258])