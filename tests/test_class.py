#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 19:59:07 2021

@author: Alexander Southan
"""

import numpy as np
import unittest

from src.pyTitration.titration import titration


class TestController(unittest.TestCase):

    def test_titration(self):

        self.titration = titration(
            [[1E7, 0, 0]], [[1E-16]], [0.1], [0.1], [1], [0])
