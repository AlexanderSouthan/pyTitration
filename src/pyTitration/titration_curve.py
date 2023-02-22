# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:37:50 2023

@author: southan

Equations are from:
https://chemistry.stackexchange.com/questions/136188/which-make-hco3-to-show-two-ph-values-at-two-scenarios/136203#136203
"""

import numpy as np
from sympy import Symbol
from sympy import solveset
import matplotlib.pyplot as plt

from little_helpers.array_tools import closest_index

class titration():
    def __init__(self, k, mode='acid', kw=1E-14):

        self.set_basic_params(k, mode=mode, kw=kw)
        self.latest_curve = None

    def set_basic_params(self, k, mode='acid', kw=1E-14):
        self.k = k
        self.mode = mode
        self.kw = kw

        v_titrant = Symbol('v_titrant')
        c_titrant = Symbol('c_titrant')
        v_analyte = Symbol('v_analyte')
        c_analyte = Symbol('c_analyte')
        c_analytix = Symbol('h_plus') if mode == 'acid' else self.kw/Symbol('h_plus')

        self.equation = (v_titrant*(c_titrant + c_analytix - kw/c_analytix) -
                         v_analyte*(c_analyte*(k[0]*c_analytix**2+2*k[0]*k[1]*c_analytix+3*k[0]*k[1]*k[2])/
                                    (c_analytix**3+k[0]*c_analytix**2+k[0]*k[1]*c_analytix+k[0]*k[1]*k[2]) - c_analytix + kw/c_analytix))
        # self.latest_curve = None

    def curve(self, v_analyte, c_analyte, c_titrant, data_points=100,
              indep_var='v_titrant', indep_var_min=0, indep_var_max=10):
        if indep_var == 'v_titrant':
            v_titrant = np.linspace(indep_var_min, indep_var_max, data_points)
            c_h_plus = np.empty_like(v_titrant)
            for idx, curr_v in enumerate(v_titrant):
                c_h_plus[idx] = max(solveset(self.equation.subs(
                    {'v_titrant': curr_v, 'v_analyte': v_analyte,
                      'c_analyte': c_analyte, 'c_titrant': c_titrant}), 'h_plus'))
            ph = -np.log10(c_h_plus)
        elif indep_var == 'pH':
            ph = np.linspace(indep_var_min, indep_var_max, data_points)
            c_h_plus = 10**(-ph)
            v_titrant = np.empty_like(ph)
            for idx, curr_c in enumerate(c_h_plus):
                v_titrant[idx] = max(solveset(self.equation.subs(
                    {'h_plus': curr_c, 'v_analyte': v_analyte,
                      'c_analyte': c_analyte, 'c_titrant': c_titrant}), 'v_titrant'))
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

if __name__ == "__main__":
    acid_titration = titration([4.46E-7, 4.69E-11, 0], mode='acid')
    x3, y3 = acid_titration.curve(0.01, 0.1, 0.1, indep_var_min=0, indep_var_max=0.03)
    
    base_titration = titration([2.13E-4, 2.24E-8, 0], mode='base')
    x4, y4 = base_titration.curve(0.01, 0.1, 0.1, indep_var_min=0, indep_var_max=0.03)
    x5, y5 = base_titration.curve(0.01, 0.1, 0.1, indep_var='pH', indep_var_min=4, indep_var_max=11)
    
    # plt.plot(x5, y5)
    plt.plot(x3, y3, x4, y4, x5, y5)
