# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 16:35:45 2023

@author: southan
"""

import numpy as np
import pandas as pd

k_values = {}
k_values['acid'] = pd.DataFrame([], index=['k_1', 'k_2', 'k_3'])

acids = []
acids.append(pd.Series([10**7, 0, 0], name='hydrochloric acid', index=k_values['acid'].index))
acids.append(pd.Series([10**-2.14, 10**-7.2, 10**-12.37], name='phosphoric acid', index=k_values['acid'].index))
acids.append(pd.Series([10**2.8, 10**-1.99, 0], name='sulfuric acid', index=k_values['acid'].index))
acids.append(pd.Series([10**-4.75, 0, 0], name='acetic acid', index=k_values['acid'].index))
acids.append(pd.Series([10**-3.13, 10**-4.76, 10**-6.39], name='citric acid', index=k_values['acid'].index))
acids.append(pd.Series([4.46E-7, 4.69E-11, 0], name='carbonic acid', index=k_values['acid'].index))
acids.append(pd.Series([10**-6.15, 0, 0], name='2-(N-morpholino)ethanesulfonic acid (MES)', index=k_values['acid'].index))
acids.append(pd.Series([10**-3, 10**-7.48, 0], name='4-(2-hydroxyethyl)-1-piperazineethanesulfonic acid (HEPES)', index=k_values['acid'].index))
acids.append(pd.Series([10**-9.25, 0, 0], name='ammonium chloride', index=k_values['acid'].index))
acids.append(pd.Series([10**-15.74, 0, 0], name='water', index=k_values['acid'].index))

k_values['acid'] = pd.concat([k_values['acid']] + acids, axis=1)

k_values['base'] = pd.DataFrame().reindex_like(k_values['acid'])

k_values['base'][k_values['acid'] > 0] = 10**-(14 + np.log10(k_values['acid'][k_values['acid'] > 0]))
k_values['base'].fillna(0, inplace=True)
k_values['base'] = k_values['base'].transform(np.sort).iloc[::-1]