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
        'idx_e': [(2, 2), (3, 3)],
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
        'y_colors': ['r', 'r', 'b', 'b', 'k', 'g', 'g', 'g'],
        'y_sizes': [2.0] * 2 + [1.0] * 5,
        'y_styles': ['-', '--'] * 2 + [':'] + ['-', '--', ':'],
        'y_name': '$n_{b}$',
        'v_label_color': 'r',
        'v_label': '$- 10 \\mathrm{log}_{10} \\langle \\tilde{Q}^{2} \\rangle$',
        'v_limits': [-2.5, 17.5],
        'v_ticks': [i * 5 for i in range(4)],
        'v_ticks_minor': [i * 1 - 3 for i in range(20)],
        'v_twin_label_color': 'g',
        'v_twin_label': '$\\langle \\tilde{\\beta}^{\\dagger} \\tilde{\\beta} \\rangle$',
        'v_twin_scale': 'log',
        'v_twin_limits': [10**-4, 10**4],
        'v_twin_tick_labels': [('$10^{' + str(i * 2 - 3) + '}$') for i in range(4)],
        'v_twin_ticks': [10.0**(i * 2 - 3) for i in range(4)],
        'v_twin_ticks_minor': [10.0**(i * 2 - 5) for i in range(20)],
        'show_legend': True,
        'legend_location': 'upper left',
        'height': 4.8,
        'width': 9.6,
        'label_font_size': 32,
        'legend_font_size': 28,
        'tick_font_size': 28
    }
}

# function to calculate the ratio and variance
def func_rat_var_n_beta(system_params, val, logger, results):
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
    # get measure dynamics
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # variance
    var = np.min(np.transpose(M)[0])
    # calculate hyperbolic angles
    rat = G_p / G_0
    r = np.arctanh(rat)
    chr = np.cosh(r)
    shr = np.sinh(r)
    # get phonon number in the Bogoluibov mode
    n_beta = np.mean([(chr**2 + shr**2) * (v[0] + v[1] - 1) / 2 + shr**2 + chr * shr * (v[0] - v[1]) for v in M])
    # update results
    results.append((val, [rat, var, n_beta]))

if __name__ == '__main__':
    # low thermal phonons with RWA
    params['system']['ns'][1] = 10.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_n_beta, looper='XLooper', file_path_prefix='data/v2.2-qom-v0.9.0/4_rwa_n=10')
    rat, var_0_rwa, n_beta_0_rwa = np.transpose(looper.results['V']).tolist()

    # high thermal phonons with RWA
    params['system']['ns'][1] = 1000.0
    params['system']['t_rwa'] = True
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_n_beta, looper='XLooper', file_path_prefix='data/v2.2-qom-v0.9.0/4_rwa_n=1000')
    _, var_1_rwa, n_beta_1_rwa = np.transpose(looper.results['V']).tolist()

    # low thermal phonons without RWA
    params['system']['ns'][1] = 10.0
    params['system']['t_rwa'] = False
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_n_beta, looper='XLooper', file_path_prefix='data/v2.2-qom-v0.9.0/4_wrwa_n=10')
    _, var_0_wrwa, _ = np.transpose(looper.results['V']).tolist()

    # high thermal phonons without RWA
    params['system']['ns'][1] = 1000.0
    params['system']['t_rwa'] = False
    looper = run_loopers_in_parallel(SystemClass=None, params=params, func=func_rat_var_n_beta, looper='XLooper', file_path_prefix='data/v2.2-qom-v0.9.0/4_wrwa_n=1000')
    _, var_1_wrwa, _ = np.transpose(looper.results['V']).tolist()

    # plotter
    plotter = MPLPlotter(axes={
        'X': rat,
        'Y': [10, 1000]
    }, params=params['plotter'])
    plotter.update(xs=rat, vs=- 10 * np.log10([var_0_rwa, var_1_rwa, var_0_wrwa, var_1_wrwa, [0.5] * len(var_0_rwa)]))
    ax = plotter.update_twin_axis(xs=rat, vs=[n_beta_0_rwa, n_beta_1_rwa, [1.0] * len(n_beta_0_rwa)])
    plotter.show(True)