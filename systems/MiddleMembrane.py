#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate Membrane-in-the-middle systems."""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.0'
__created__ = '2022-04-01'
__updated__ = '2023-06-26'

# dependencies
import numpy as np
import scipy.integrate as si

# qom modules
from qom.systems import BaseSystem

class MM_01(BaseSystem):
    r"""Class to simulate a membrane-in-the-middle system driven by a modulated laser using constant mode amplitudes.

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        alphas      (list) base and sideband amplitudes of the optical mode :math:`[ \alpha_{0}, \alpha_{-}, \alpha_{+} ]`. Default is :math:`[ 2.0, 0.8, 0.8 ]`.
        betas       (list) base and sideband amplitudes of the mechanical mode :math:`[ \beta_{0}, \beta_{-}, \beta_{+} ]`. Default is :math:`[ 100.0, 25.0, 62.5 ]`.
        Delta_norm  (float) normalized effective detuning of the cavity from the laser :math:`\Delta / \omega_{m}`. Default is :math:`1.0`.
        g_norm      (float) normalized optomechanical coupling strength :math:`g / \omega_{m}`. Default is :math:`10^{-4}`.
        gamma_norm  (float) normalized mechanical damping rate :math:`\gamma / \omega_{m}`. Default is :math:`10^{-6}`.
        kappa_norm  (float) normalized optical decay rate :math:`\kappa / \omega_{m}`. Default is :math:`0.1`.
        ns          (list) quanta of thermal photons and phonons :math:`[ n_{a}, n_{b} ]`. Default is :math:`[ 0.0, 10.0 ]`.
        Omega_norms (list) normalized modulation frequencies :math:`[ \Omega_{a} / \omega_{m}, \Omega_{b} / \omega_{m} ]`. Default is :math:`[ 2.0, 2.0 ]`.
        t_rwa       (bool) option to work under RWA. Default is `True`.
        ========    ============================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default looping parameters
    looper_defaults = {
        'alpha_0': {
            'var': 'alphas',
            'idx': 0,
            'min': 0.0,
            'max': 4.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'alpha_m': {
            'var': 'alphas',
            'idx': 1,
            'min': 0.0,
            'max': 2.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'alpha_p': {
            'var': 'alphas',
            'idx': 2,
            'min': 0.0,
            'max': 2.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'beta_0': {
            'var': 'betas',
            'idx': 0,
            'min': 0.0,
            'max': 200.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'beta_m': {
            'var': 'betas',
            'idx': 1,
            'min': 0.0,
            'max': 100.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'beta_p': {
            'var': 'betas',
            'idx': 2,
            'min': 0.0,
            'max': 100.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'Delta_norm': {
            'var': 'Delta_norm',
            'min': -2.0,
            'max': 2.0,
            'dim': 1001,
            'scale': 'linear'
        },
        'g_norm': {
            'var': 'g_norm',
            'min': 1e-6,
            'max': 1e-2,
            'dim': 1001,
            'scale': 'log'
        },
        'gamma_norm': {
            'var': 'gamma_norm',
            'min': 1e-8,
            'max': 1e-4,
            'dim': 1001,
            'scale': 'log'
        },
        'kappa_norm': {
            'var': 'kappa_norm',
            'min': 1e-2,
            'max': 1e2,
            'dim': 1001,
            'scale': 'log'
        },
        'n_a': {
            'var': 'ns',
            'idx': 0,
            'min': 1e-2,
            'max': 1e2,
            'dim': 1001,
            'scale': 'log'
        },
        'n_b': {
            'var': 'ns',
            'idx': 1,
            'min': 1e-1,
            'max': 1e3,
            'dim': 1001,
            'scale': 'log'
        },
        'Omega_a_norm': {
            'var': 'Omega_norms',
            'idx': 0,
            'min': 1.5,
            'max': 2.5,
            'dim': 1001,
            'scale': 'linear'
        },
        'Omega_b_norm': {
            'var': 'Omega_norms',
            'idx': 1,
            'min': 1.5,
            'max': 2.5,
            'dim': 1001,
            'scale': 'linear'
        },
        't_rwa': {
            'var': 't_rwa',
            'val': [True, False]
        }
    }

    # default system parameters
    system_defaults = {
        'alphas'        : [2.0, 0.8, 0.8],
        'betas'         : [100.0, 25.0, 62.5],
        'Delta_norm'    : 1.0,
        'g_norm'        : 1e-4,
        'gamma_norm'    : 1e-6,
        'kappa_norm'    : 0.1,
        'ns'            : [0.0, 10.0],
        'Omega_norms'   : [2.0, 2.0],
        't_rwa'         : True
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for MM_01."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='MM_01',
            desc='Modulated Membrane-in-the-middle System',
            num_modes=2,
            cb_update=cb_update
        )

        # set drift matrix as constant under RWA
        assert type(self.params['t_rwa']) is bool, 'Parameter "t_rwa" should be of type boolean'
        self.is_A_constant = self.params['t_rwa']

    def get_A(self, modes, c, t):
        """Method to obtain the drift matrix.

        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        c : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the values are calculated.
        
        Returns
        -------
        A : numpy.ndarray
            Drift matrix.
        """

        # extract frequently used variables
        g_norm                      = self.params['g_norm']
        gamma_norm                  = self.params['gamma_norm']
        kappa_norm                  = self.params['kappa_norm']
        Omega_a_norm, Omega_b_norm  = self.params['Omega_norms']

        # with RWA
        if self.params['t_rwa']:
            # normalized effective couplings
            G_minus_norm, G_plus_norm, G_tilde_minus_norm, G_tilde_plus_norm = self.get_params_G_norms(
                c=c
            )

            # X quadratures
            self.A[0][0]    = - kappa_norm / 2.0
            self.A[0][3]    = - G_minus_norm
            # Y quadratures
            self.A[1][1]    = - kappa_norm / 2.0
            self.A[1][2]    = G_plus_norm
            # Q quadratures
            self.A[2][1]    = - G_minus_norm
            self.A[2][2]    = - gamma_norm / 2.0
            self.A[2][3]    = - G_tilde_minus_norm
            # P quadratures
            self.A[3][0]    = G_plus_norm
            self.A[3][2]    = G_tilde_plus_norm
            self.A[3][3]    = - gamma_norm / 2.0

        # without RWA
        else:
            # extract frequently used variables
            alpha_0, alpha_m, alpha_p   = self.params['alphas']
            beta_0, beta_m, beta_p      = self.params['betas']
            Delta_norm                  = self.params['Delta_norm']
            
            # modes
            alpha   = alpha_0 + alpha_m * np.exp(1.0j * Omega_a_norm * t) + alpha_p * np.exp(-1.0j * Omega_a_norm * t)
            beta    = beta_0 + beta_m * np.exp(1.0j * Omega_b_norm * t) + beta_p * np.exp(-1.0j * Omega_b_norm * t)

            # X quadratures
            self.A[0][0]    = - kappa_norm / 2.0
            self.A[0][1]    = Delta_norm
            self.A[0][2]    = - 8.0 * g_norm * np.real(beta) * np.imag(alpha)
            # Y quadratures
            self.A[1][0]    = - Delta_norm
            self.A[1][1]    = - kappa_norm / 2.0
            self.A[1][2]    = 8.0 * g_norm * np.real(beta) * np.real(alpha)
            # Q quadratures
            self.A[2][2]    = - gamma_norm / 2.0
            self.A[2][3]    = 1.0
            # P quadratures
            self.A[3][0]    = 8.0 * g_norm * np.real(beta) * np.real(alpha)
            self.A[3][1]    = 8.0 * g_norm * np.real(beta) * np.imag(alpha)
            self.A[3][2]    = - 1.0 + 4.0 * g_norm * np.real(np.conjugate(alpha) * alpha)
            self.A[3][3]    = - gamma_norm / 2.0

        return self.A

    def get_params_ratio(self, c):
        """Method to obtain the squeezing ratio.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        ratio : numpy.ndarray
            Squeezing ratio.
        """

        # get normalized effective couplings
        G_minus_norm, G_plus_norm, _, _ = self.get_params_G_norms(
            c=c
        )

        return (G_plus_norm - G_minus_norm) / (G_plus_norm + G_minus_norm)

    def get_params_G_norms(self, c):
        """Method to obtain the normalized effective couplings.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        values : numpy.ndarray
            Normalized effective couplings in the order G_minus, G_plus, G_tilde_plus, G_tilde_minus.
        """

        # extract frequently used parameters
        alpha_0, alpha_m, alpha_p   = self.params['alphas']
        beta_0, beta_m, beta_p      = self.params['betas']
        g_norm                      = self.params['g_norm']
        
        # normalized effective couplings
        G_0_norm        = 2 * g_norm * (2 * alpha_0 * beta_0 + (alpha_m + alpha_p) * (beta_m + beta_p))
        G_1_norm        = 2 * g_norm * (alpha_0 * (beta_m + beta_p) + 2 * alpha_p * beta_0)
        G_tilde_0_norm  = 2 * g_norm * (alpha_0**2 + alpha_m**2 + alpha_p**2)
        G_tilde_1_norm  = 2 * g_norm * (alpha_0 * (alpha_m + alpha_p))

        # substituted expressions
        G_minus_norm        = G_0_norm - G_1_norm
        G_plus_norm         = G_0_norm + G_1_norm
        G_tilde_minus_norm  = G_tilde_0_norm - G_tilde_1_norm
        G_tilde_plus_norm   = G_tilde_0_norm + G_tilde_1_norm

        return G_minus_norm, G_plus_norm, G_tilde_minus_norm, G_tilde_plus_norm

    def get_D(self, modes, corrs, c, t):
        """Method to obtain the noise matrix.
        
        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        corrs : numpy.ndarray
            Quantum correlations.
        c : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the values are calculated.
        
        Returns
        -------
        D : numpy.ndarray
            Noise matrix.
        """

        # extract frequently used variables
        gamma_norm  = self.params['gamma_norm']
        kappa_norm  = self.params['kappa_norm']
        n_a, n_b    = self.params['ns']
        
        # optical mode
        self.D[0][0]    = kappa_norm * (n_a + 0.5)
        self.D[1][1]    = kappa_norm * (n_a + 0.5)
        # mechanical mode
        self.D[2][2]    = gamma_norm * (n_b + 0.5)
        self.D[3][3]    = gamma_norm * (n_b + 0.5)

        return self.D
    
    def get_ivc(self):
        """Method to obtain the initial values of the modes, correlations and derived constants and controls.
        
        Returns
        -------
        iv_modes : numpy.ndarray
            Initial values of the classical modes.
        iv_corrs : numpy.ndarray
            Initial values of the quantum correlations.
        c : numpy.ndarray
            Derived constants and controls.
        """

        # extract frequently used variables
        n_a, n_b = self.params['ns']
 
        # initial mode values
        iv_modes = np.zeros(self.num_modes, dtype=np.complex_)

        # initial quadrature correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = n_a + 0.5
        iv_corrs[1][1]  = n_a + 0.5
        iv_corrs[2][2]  = n_b + 0.5
        iv_corrs[3][3]  = n_b + 0.5

        return iv_modes, iv_corrs, np.empty(0)

    def get_mode_rates(self, modes, c, t):
        """Method to obtain the rates of change of the modes.

        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        c : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the values are calculated.
        
        Returns
        -------
        mode_rates : numpy.ndarray
            Rate of change of the modes.
        """

        return np.array([0.0, 0.0], dtype=np.complex_)

    def get_var_Q_ss_rwa(self, c):
        """Method to obtain the steady-state variance of the position quadrature under RWA.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        var_Q_ss_rwa : float
            Variance of the position quadrature.
        """

        # extract frequently used variables
        gamma_norm  = self.params['gamma_norm']
        kappa_norm  = self.params['kappa_norm']
        n_a, n_b    = self.params['ns']

        # get normalized effective couplings
        G_minus_norm, G_plus_norm, G_tilde_minus_norm, _ = self.get_params_G_norms(
            c=c
        )

        # get squeezing ratio
        ratio = self.get_params_ratio(
            c=c
        )

        # substituted expressions
        r = np.arctanh(ratio)
        h = 2.0 * G_plus_norm * G_minus_norm / kappa_norm + gamma_norm / 2.0

        # corner cases
        if G_minus_norm * h == 0.0:
            _coeff = np.inf
        else:
            _coeff = G_tilde_minus_norm * G_plus_norm / (G_minus_norm * h)
        if G_tilde_minus_norm**2 - h**2 == 0.0:
            return np.inf
            
        # steady-state variance
        return h * np.exp(- 2.0 * r) / 2.0 / (G_tilde_minus_norm**2 - h**2) * (gamma_norm * (n_b + 0.5) * (_coeff * np.exp(-2.0 * r) - np.exp(2.0 * r)) - 4.0 * G_plus_norm * G_minus_norm / kappa_norm * (n_a + 0.5) * (1.0 + _coeff))

    def get_var_Q_ft_rwa(self, c):
        """Method to obtain the variance of the position quadrature using the Fourier transform under RWA.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        var_Q_ft_rwa : float
            Variance of the position quadrature.
        """

        # extract frequently used variables
        gamma_norm  = self.params['gamma_norm']
        kappa_norm  = self.params['kappa_norm']
        n_a, n_b    = self.params['ns']

        # normalized effective couplings
        G_minus_norm, G_plus_norm, G_tilde_minus_norm, G_tilde_plus_norm = self.get_params_G_norms(
            c=c
        )

        # expressions
        _den    = lambda omega_norm: (kappa_norm - 2.0j * omega_norm)**2 * (- 4.0 * G_tilde_minus_norm * G_tilde_plus_norm + (2.0 * omega_norm + 1.0j * gamma_norm)**2) + 8.0 * G_minus_norm * G_plus_norm * (2.0 * omega_norm + 1.0j * kappa_norm) * (2.0 * omega_norm + 1.0j * gamma_norm) - 16.0 * G_minus_norm**2 * G_plus_norm**2
        A_num   = lambda omega_norm: 8.0 * G_plus_norm * np.sqrt(kappa_norm) * G_tilde_minus_norm * (kappa_norm - 2.0j * omega_norm)
        B_num   = lambda omega_norm: 4.0 * G_minus_norm * np.sqrt(kappa_norm) * (- 4.0 * G_minus_norm * G_plus_norm + (2.0 * omega_norm + 1.0j * kappa_norm) * (2.0 * omega_norm + 1.0j * gamma_norm))
        C_num   = lambda omega_norm: 2.0 * (kappa_norm - 2.0j * omega_norm) * np.sqrt(gamma_norm) * (4.0 * G_minus_norm * G_plus_norm + (kappa_norm - 2.0j * omega_norm) * (gamma_norm - 2.0j * omega_norm))
        D_num   = lambda omega: - 4.0 * G_tilde_minus_norm * (kappa_norm - 2.0j * omega)**2 * np.sqrt(gamma_norm)

        # fluctuation spectrum
        S_Q = lambda omega_norm: (A_num(- omega_norm) * A_num(omega_norm) + B_num(- omega_norm) * B_num(omega_norm)) / _den(- omega_norm) / _den(omega_norm) * (n_a + 0.5) + (C_num(- omega_norm) * C_num(omega_norm) + D_num(- omega_norm) * D_num(omega_norm)) / _den(- omega_norm) / _den(omega_norm) * (n_b + 0.5)

        # variance
        return 1.0 / 2.0 / np.pi * si.quad(S_Q, -np.inf, np.inf)[0]