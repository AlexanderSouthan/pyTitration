# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:37:50 2023

@author: southan
"""

import numpy as np
from sympy import Symbol
from sympy import solveset
from scipy.optimize import brentq

from little_helpers.array_tools import closest_index
from little_helpers.num_derive import derivative


class solution():
    def __init__(self, k_solutes, c_solutes, prot_left, kw=1E-14):
        k_len = [len(curr_k_set) for curr_k_set in k_solutes]
        if np.ptp(k_len) != 0:
            max_len = max(k_len)
            k_solutes = [curr_k_set + [0]*(max_len-len(curr_k_set))
                         for curr_k_set in k_solutes]
        self.k_solutes = np.asarray(k_solutes)
        self.c_solutes = np.asarray(c_solutes)
        self.prot_left = np.asarray(prot_left, dtype=int)
        self.kw = kw

        self.h_plus = Symbol('h_plus')

        self._calc_equation()

    def calc_ph(self, ph_min=0, ph_max=16):
        c_h_plus = brentq(
            self._calc_equation_value, 10**-ph_max, 10**-ph_min, xtol=1E-16)
        return -np.log10(c_h_plus)

    def _calc_f(self):
        n = np.sum(self.k_solutes!=0, axis=1)
        assert self.k_solutes.shape[0] == len(self.prot_left), (
            'Number of dissociation constant sets and entries in prot_left '
            'must be equal, but are {} and {}.').format(
                self.k_solutes.shape[0], len(self.prot_left))
        assert (self.prot_left <= n).all(), (
            'prot_left is {} and n is {}, but prot_left must be smaller '
            'or equal n.').format(self.prot_left, n)

        funcs = []
        for curr_prot, curr_n, curr_k_set in zip(self.prot_left, n,
                                                 self.k_solutes):
            m_range = np.arange(0, curr_n+1)
            mask = np.ones_like(m_range, dtype=bool)
            mask[curr_prot] = False

            m_range = m_range[mask]
            front_factors = m_range - curr_prot

            alpha = []
            for m in m_range:
                numer = self.h_plus**m * np.prod(curr_k_set[:curr_n-m])
                denom = self.h_plus**curr_n
                for ii, _ in enumerate(curr_k_set):
                    denom += self.h_plus**(curr_n-ii-1)*np.prod(
                        curr_k_set[:ii+1])
                alpha.append(numer/denom)

            func = 0
            for curr_alpha, curr_ff in zip(alpha, front_factors):
                func += curr_ff * curr_alpha
            funcs.append(func)
        return np.asarray(funcs)

    def _calc_equation(self):
        delta = self.h_plus - self.kw/self.h_plus
        func = self._calc_f()
        self.equation = np.sum(func*self.c_solutes, axis=0)+delta
        return self.equation

    def _calc_equation_value(self, c_h_plus):
        return float(self.equation.subs({'h_plus': c_h_plus}))


class titration():
    def __init__(self, k_analyte, k_titrant, c_analyte, c_titrant,
                 prot_left_ana, prot_left_tit, kw=1E-14):
        """
        Initialize a titration instance.

        Parameters
        ----------
        k_analyte : list or ndarray
            Contains the acid dissociation constants (K_a) of the solutes in
            the solution to be analyzed, must be in decreasing order for each
            solute. The form is aa list containing lists with the K_a values.
            For example, the monoprotic acid hydrochloric acid would be
            [[1E7]] and aa mixture of hydrochloric acid and phosphoric acid
            [[1E7], [7.24E-3, 6.31E-8, 4.27E-13]]. The numbers must be the acid
            dissociation constants also if the solutes are in their basic form!
        k_titrant : list or ndarray
            Contains the acid dissociation constanst of the solutes in the
            titrant, i.e. the solution added to the analyte during titration.
            The requirements for this list are the same as for k_analyte.
        c_analyte : ndarray or list of float
            A list containing the concentrations in mol/L of the solutes in the
            analyte solution. For each set of K_a values provided in k_analyte,
            one value in c_analyte is needed.
        c_titrant : TYPE
            A list containing the concentrations in mol/L of the solutes in the
            titrant solution. For each set of K_a values provided in k_titrant,
            one value in c_titrant is needed.
        prot_left_ana : ndarray or list of float
            A list containing the number of residual acidic protons (for which
            K_a values are provided in k_analyte) of each solute in the analyte
            solution. For example, H3PO4 would be 3, H2PO4- would be 2, HCl
            would be 1, NO3- would be 0, and OH- would be 0.
        prot_left_tit : ndarray or list of float
            A list containing the number of residual acidic protons (for which
            K_a values are provided in k_titrant) of each solute in the titrant
            solution.
        kw : float, optional
            The ion product of pure water. The default is 1E-14.

        Returns
        -------
        None.

        """

        self.set_basic_params(k_analyte, k_titrant, c_analyte, c_titrant,
                              prot_left_ana, prot_left_tit, kw=kw)
        self.latest_curve = None

    def set_basic_params(self, k_analyte, k_titrant, c_analyte, c_titrant,
                         prot_left_ana, prot_left_tit, kw=1E-14):
        self.analyte = solution(k_analyte, c_analyte, prot_left_ana, kw=kw)
        self.titrant = solution(k_titrant, c_titrant, prot_left_tit, kw=kw)
        self.ph_bounds = np.sort(
            [self.analyte.calc_ph(), self.titrant.calc_ph()])
        self.h_plus_bounds = (10**-self.ph_bounds)[::-1]

        self._calc_equation()
        # self.latest_curve = None

    def curve(self, v_analyte, data_points=100,
              indep_var='pH', indep_var_min=0, indep_var_max=10):
        """
        Calulate a titration curve of the mixture given by the arguments for
        the init method.

        Calculation of titration curves is done according to:
        Analytical Chemistry 1996, 68 (4), 585-590. DOI: 10.1021/ac950430l.
        This allows to calculate titration curves of arbitrary mixtures.

        Parameters
        ----------
        v_analyte : float
            The volume of the analyte solution in litres.
        data_points : floaat, optional
            The number of data points of the calculated titration curve. The
            default is 100.
        indep_var : str, optional
            Can either be 'v_titrant' or 'pH'. If it is 'v_titrant', the
            titration curve equation has to be solved for the H+ concentration
            which is not an easy task. This does work in principle at the
            moment, but the bounds used could be improved. The default is 'pH'
            and the equation is solved for the titrant volume which is very
            easy.
        indep_var_min : float, optional
            The minimum value of the independent variable used for the
            calculation, as defined by indep_var. The default is 0.
        indep_var_max : TYPE, optional
            The maximum value of the independent variable used for the
            calculation, as defined by indep_var.The default is 10.

        Returns
        -------
        tuple of ndarrays
            The titration curve. The first element is the volume of the
            titrant, the second element the resulting pH.

        """
        if indep_var == 'v_titrant':
            v_titrant = np.linspace(indep_var_min, indep_var_max, data_points)
            c_h_plus = np.empty_like(v_titrant)
            for idx, curr_v in enumerate(v_titrant):
                c_h_plus[idx] = brentq(
                    self._calc_equation_value, 0,  # self.h_plus_bounds[0],
                    0.99*self.h_plus_bounds[1], args=(v_analyte, curr_v),
                    xtol=1E-16)
            ph = -np.log10(c_h_plus)
        elif indep_var == 'pH':
            if (indep_var_min <= self.ph_bounds[0]) or (
                    indep_var_max >= self.ph_bounds[1]):
                raise ValueError(
                    'The minimum and maximum pH values must be between {} and '
                    '{} (the analyte/titrant pH values), but are {} '
                    'and {}.'.format(self.ph_bounds[0], self.ph_bounds[1],
                                     indep_var_min, indep_var_max))

            ph = np.linspace(indep_var_min, indep_var_max, data_points)
            c_h_plus = 10**(-ph)
            v_titrant = np.empty_like(ph)

            sol = solveset(self.equation.subs(
                {'v_analyte': v_analyte}), 'v_titrant')
            for idx, curr_c in enumerate(c_h_plus):
                v_titrant[idx] = max(sol.subs({'h_plus': curr_c}))
        else:
            raise ValueError('indep_var must either be \'v_titrant\' or'
                             ' \'pH\', but is \'{}\'.'.format(indep_var))

        self.latest_curve = (v_titrant, ph)
        return self.latest_curve

    def curve_derivative(self, order=1):
        if self.latest_curve is None:
            raise ValueError('self.latest_curve is None. Run self.curve(...) '
                             'first or provide titration curve data manually.')

        deriv = derivative(self.latest_curve[0], [self.latest_curve[1]],
                           order=order)
        return (self.latest_curve[0], deriv[0])

    def volume_between_ph(self, ph1, ph2):
        if self.latest_curve is None:
            raise ValueError('self.latest_curve is None. Run self.curve(...) '
                             'first or provide titration curve data manually.')

        idx = closest_index([ph1, ph2], self.latest_curve[1])
        volume = np.diff(self.latest_curve[0][idx]).item()
        return volume

    def export_titration_curve(self, file_name):
        if self.latest_curve is None:
            raise ValueError(
                'No titration curve available for export, calculate one '
                'first.')
        else:
            exp = np.asarray(self.latest_curve).T
            np.savetxt(file_name + '.csv', exp, delimiter=',')

    def _calc_equation(self):
        v_titrant = Symbol('v_titrant')
        v_analyte = Symbol('v_analyte')

        self.equation = -v_titrant/v_analyte - (
            self.analyte.equation/self.titrant.equation)

    def _calc_equation_value(self, c_h_plus, v_analyte, v_titrant):
        return float(self.equation.subs(
            {'h_plus': c_h_plus, 'v_analyte': v_analyte,
             'v_titrant': v_titrant}))
