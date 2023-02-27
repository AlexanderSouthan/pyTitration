# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:27:17 2023

@author: southan
"""

import matplotlib.pyplot as plt

from pyTitration import titration
from pyTitration import k_values

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
    0.5, indep_var='v_titrant', indep_var_min=0, indep_var_max=1)
x1, y1 = acid_titration.curve(
    0.5, indep_var_min=3.68, indep_var_max=12, indep_var='pH')
x2, y2 = acid_mix_titration.curve(
    0.5, indep_var='pH', indep_var_min=2, indep_var_max=12.8)
x3, y3 = base_titration.curve(
    0.5, indep_var='pH', indep_var_min=2, indep_var_max=11.65)

plt.plot(x0, y0, x1, y1, x2, y2, x3, y3)
plt.xlabel('$V_\mathrm{titrant}$ [L]')
plt.ylabel('pH')