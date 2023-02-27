# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:37:50 2023

@author: southan
"""

import numpy as np
from sympy import Symbol
from sympy import solveset

from little_helpers.array_tools import closest_index

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
        self.k_analyte = np.asarray(k_analyte)
        self.k_titrant = np.asarray(k_titrant)
        self.c_analyte = np.asarray(c_analyte)
        self.c_titrant = np.asarray(c_titrant)
        self.prot_left_ana = prot_left_ana
        self.prot_left_tit = prot_left_tit
        self.kw = kw

        self.h_plus = Symbol('h_plus')

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
            which is not an easy task. This basically does not work aat the
            moment. The default is 'pH' and the equation is solved for the
            titrant volume which is very easy.
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
            # this part is basically unfunctional at the moment, the root
            # finding procedure would have to be changed so that a solution
            # is found reliably.
            v_titrant = np.linspace(indep_var_min, indep_var_max, data_points)
            c_h_plus = np.empty_like(v_titrant)
            for idx, curr_v in enumerate(v_titrant):
                c_h_plus[idx] = max(solveset(self.equation.subs(
                    {'v_titrant': curr_v, 'v_analyte': v_analyte}), 'h_plus'))
            ph = -np.log10(c_h_plus)
        elif indep_var == 'pH':
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

    def volume_between_ph(self, ph1, ph2):
        if self.latest_curve is None:
            raise ValueError('self.latest_curve is None. Run self.curve(...) '
                             'first or provide titration curve data manually.')

        idx = closest_index(sorted([ph1, ph2]), self.latest_curve[1])
        volume = np.diff(self.latest_curve[0][idx]).item()
        return volume

    def _calc_equation(self):
        v_titrant = Symbol('v_titrant')
        v_analyte = Symbol('v_analyte')

        delta = self.h_plus - self.kw/self.h_plus

        analyte_func = self._calc_f(
            self.k_analyte, self.prot_left_ana)
        titrant_func = self._calc_f(
            self.k_titrant, self.prot_left_tit)

        self.equation = -v_titrant/v_analyte - (
            np.sum(analyte_func*self.c_analyte, axis=0)+delta)/(
                np.sum(titrant_func*self.c_titrant, axis=0)+delta)

    def _calc_f(self, k, protons_left):
        n = np.sum(k!=0, axis=1)
        assert k.shape[0] == len(protons_left), (
            'Number of dissociation constant sets and entries in protons_left '
            'must be equal, but are {} and {}.').format(
                k.shape[0], len(protons_left))
        assert (protons_left <= n).all(), (
            'protons_left is {} and n is {}, but protons_left must be smaller '
            'or equal n.').format(protons_left, n)

        funcs = []
        for curr_prot, curr_n, curr_k_set in zip(protons_left, n, k):
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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from k_values import k_values

    simple_titration = titration(
        k_analyte=[k_values['acid']['water']],
        k_titrant=[k_values['acid']['hydrochloric acid']],
        c_analyte=[0.1], c_titrant=[0.1],
        prot_left_ana=[1], prot_left_tit=[1])
    acid_titration = titration(
        k_analyte=[k_values['acid']['carbonic acid']],
        k_titrant=[k_values['acid']['water']],
        c_analyte=[0.1], c_titrant=[0.1],
        prot_left_ana=[2], prot_left_tit=[0])
    base_titration = titration(
        k_analyte=[k_values['acid']['carbonic acid']],
        k_titrant=[k_values['acid']['hydrochloric acid']],
        c_analyte=[0.1], c_titrant=[0.1],
        prot_left_ana=[0], prot_left_tit=[1])
    acid_mix_titration = titration(
        k_analyte=[[10**-2.87], [10**-4.8], [10**-6.96]],
        k_titrant=[k_values['acid']['water']],
        c_analyte=[0.1, 0.1, 0.1],
        c_titrant=[0.3], prot_left_ana=[1, 1, 1],
        prot_left_tit=[0])

    # x0, y0 = simple_titration.curve(
    #     0.5, indep_var='v_titrant', indep_var_min=0, indep_var_max=1)
    x1, y1 = acid_titration.curve(
        0.5, indep_var_min=2, indep_var_max=12, indep_var='pH')
    x2, y2 = acid_mix_titration.curve(
        0.5, indep_var='pH', indep_var_min=2, indep_var_max=12.8)
    x3, y3 = base_titration.curve(
        0.5, indep_var='pH', indep_var_min=2, indep_var_max=12)
    
    plt.plot(x1, y1, x2, y2, x3, y3)
