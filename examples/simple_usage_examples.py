# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:27:17 2023

@author: southan
"""

import numpy as np
import pandas as pd
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

gelatin_type_a = pd.DataFrame(np.array(
    [[0.2861456 , 0.15787344, 0.14439366, 0.32828285, 0.46925584, 0.24440408,
     3.31573114, 1.03711968, 0, 0.20658, 0.06031731, 0.09376973, 0.2226078,
     0.03311476, 0.1259166, 0.04833901, 0.25925698, 0.48335136, 1.17953515,
     0, 0.901399, 0.08077126],
     [4.05, np.nan, np.nan, np.nan, 4.45, np.nan, np.nan, np.nan, 9, np.nan,
      np.nan, np.nan, np.nan, 10, np.nan, 5.98, 10, 12, np.nan, np.nan, np.nan,
      10]]).T,
    index=['D', 'N', 'T', 'S', 'E', 'Q', 'G', 'A', 'C', 'V', 'M', 'I', 'L',
           'Y', 'F', 'H', 'K', 'R', 'P', 'W', 'Hyp', 'Hyl'],
    columns=['mmol/g', 'pka_bjellqvist']).dropna()
gelatin_type_a['c_10% [mol/L]'] = gelatin_type_a['mmol/g']*100/1000


gelatin_titration = titration(
    k_analyte=[[curr_k] for curr_k in gelatin_type_a['pka_bjellqvist']],
    k_titrant=[k_values['acid']['water']],
    c_analyte=gelatin_type_a['c_10% [mol/L]'].tolist(), c_titrant=[4],
    prot_left_ana=[1]*len(gelatin_type_a['c_10% [mol/L]']),
    prot_left_tit=[0])

# titration of a 10 % gelatin solution in 10mM phosphate buffer
gelatin_titration_buffered = titration(
    k_analyte=[[curr_k] for curr_k in gelatin_type_a['pka_bjellqvist']] + [k_values['acid']['phosphoric acid'].tolist()],
    k_titrant=[k_values['acid']['water']],
    c_analyte=gelatin_type_a['c_10% [mol/L]'].tolist() + [0.01], c_titrant=[4],
    prot_left_ana=[1]*len(gelatin_type_a['c_10% [mol/L]']) + [3],
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
x4, y4 = gelatin_titration.curve(
    0.25, indep_var='pH', indep_var_min=1.4, indep_var_max=12, data_points=1000)
x5, y5 = gelatin_titration_buffered.curve(
    0.25, indep_var='pH', indep_var_min=1.4, indep_var_max=12, data_points=1000)

simple_titration.export_titration_curve('simple_titration')
acid_titration.export_titration_curve('acid_titration')
base_titration.export_titration_curve('base_titration')
acid_mix_titration.export_titration_curve('acid_mix_titration')
gelatin_titration.export_titration_curve('gelatin_titration')
gelatin_titration_buffered.export_titration_curve('gelatin_titration_buffered')

plt.plot(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5)
plt.xlabel('$V_\mathrm{titrant}$ [L]')
plt.ylabel('pH')
