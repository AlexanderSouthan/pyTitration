#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 19:59:07 2021

@author: Alexander Southan
"""

import matplotlib.pyplot as plt
import unittest

from src.pyTitration.titration import titration
from src.pyTitration.k_values import k_values


class TestController(unittest.TestCase):

    def test_titration(self):

        simple_titration = titration(
            k_analyte=[k_values['acid']['water']],
            k_titrant=[k_values['acid']['hydrochloric acid']],
            c_analyte=[0.1], c_titrant=[0.1],
            prot_left_ana=[0], prot_left_tit=[1])
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
        
        x0, y0 = simple_titration.curve(
            0.5, indep_var='v_titrant', indep_var_min=0, indep_var_max=1,
            data_points=1000)
        x1, y1 = acid_titration.curve(
            0.5, indep_var_min=3.68, indep_var_max=12, indep_var='pH',
            data_points=1000)
        x2, y2 = acid_mix_titration.curve(
            0.5, indep_var='pH', indep_var_min=2, indep_var_max=12.8,
            data_points=1000)
        x3, y3 = base_titration.curve(
            0.5, indep_var='pH', indep_var_min=2, indep_var_max=11.65,
            data_points=1000)
        
        simple_titration.export_titration_curve('simple_titration')
        acid_titration.export_titration_curve('acid_titration')
        base_titration.export_titration_curve('base_titration')
        acid_mix_titration.export_titration_curve('acid_mix_titration')
        
        plt.plot(x0, y0, x1, y1, x2, y2, x3, y3)
        plt.xlabel('$V_\mathrm{titrant}$ [L]')
        plt.ylabel('pH')
