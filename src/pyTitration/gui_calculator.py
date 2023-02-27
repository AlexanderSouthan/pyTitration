# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:04:58 2023

@author: southan
"""

import numpy as np
import pandas as pd
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtGui import QDoubleValidator
import pyqtgraph as pg

from pyTitration import titration
from pyTitration import k_values
from little_helpers.array_tools import y_at_x

class titration_window(QtWidgets.QWidget):
    def __init__(self, set_point):
        super().__init__()

        self.set_point = set_point
        self.error = None

        self.params = pd.Series(
            ['acid', 1E7, 0, 0, 1E-14, 0.1, 0.1, 0.5, 3, 11],
            index=['acid_base', '<i>K</i><sub>1</sub>', '<i>K</i><sub>2</sub>',
                   '<i>K</i><sub>3</sub>',
                   '<i>K</i><sub>w</sub> [mol<sup>2</sup>/L<sup>2</sup>]',
                   '<i>c</i><sub>titrant</sub> [mol/L]', '<i>c</i><sub>solution</sub> [mol/L]',
                   '<i>V</i><sub>solution</sub> [L]', 'pH<sub>min</sub>',
                   'pH<sub>max</sub>'])

        self.titration = titration([self.params.iloc[1:4].tolist()],
                                   [[1E-16]], [self.params.iloc[6]],
                                   [self.params.iloc[5]], [1], [0],
                                   kw=self.params.iloc[4])

        self.init_window()
        self.define_widgets()
        self.position_widgets()
        self.connect_event_handlers()

    def init_window(self):
        self.setGeometry(500, 500, 500, 600)  # xPos,yPos,width, heigth
        self.center()  # center function is defined below
        self.setWindowTitle('Titration curve calculator')

        self.grid_container = QtWidgets.QGridLayout()
        self.setLayout(self.grid_container)

    def define_widgets(self):
        self.error_graph = pg.PlotWidget()
        self.error_graph.setBackground('w')
        self.error_graph.setLabel('left', 'pH')
        self.error_graph.setLabel('bottom', 'titrant volume [L]')
        self.error_graph.setTitle('titration curve')
        self.error_pen = pg.mkPen(color=(0, 0, 255))
        self.error_line = self.error_graph.plot(pen=self.error_pen)

        self.acid_combo = QtWidgets.QComboBox()
        self.acid_combo.addItems(
            ['acid', 'base'])

        self.acid_presets_combo = QtWidgets.QComboBox()
        self.acid_presets_combo.addItem('')
        self.acid_presets_combo.addItems(k_values['acid'].columns)

        self.labels = {}
        self.les = {}
        onlyFloat = QDoubleValidator()
        for curr_param, curr_value in self.params.items():
            if curr_param != 'acid_base':
                self.labels[curr_param] = QtWidgets.QLabel(curr_param)
                self.les[curr_param] = QtWidgets.QLineEdit(str(curr_value))
                self.les[curr_param].setValidator(onlyFloat)

        self.titcurve_btn = QtWidgets.QPushButton('Calculate titration curve')
        self.export_btn = QtWidgets.QPushButton('Export titration curve')

    def position_widgets(self):
        vlay = QtWidgets.QVBoxLayout()
        vlay.addWidget(self.error_graph)

        hlay_k = QtWidgets.QHBoxLayout()
        hlay_k.addWidget(self.acid_combo)
        for curr_param in self.params.index[1:4]:
            hlay_k.addWidget(self.labels[curr_param])
            hlay_k.addWidget(self.les[curr_param])
        vlay.addLayout(hlay_k)

        vlay.addWidget(self.acid_presets_combo)

        vlay_labels1 = QtWidgets.QVBoxLayout()
        vlay_les1 = QtWidgets.QVBoxLayout()
        for curr_param in self.params.index[4:7]:
            vlay_labels1.addWidget(self.labels[curr_param])
            vlay_les1.addWidget(self.les[curr_param])
        vlay_labels2 = QtWidgets.QVBoxLayout()
        vlay_les2 = QtWidgets.QVBoxLayout()
        for curr_param in self.params.index[7:]:
            vlay_labels2.addWidget(self.labels[curr_param])
            vlay_les2.addWidget(self.les[curr_param])
        hlay_params = QtWidgets.QHBoxLayout()
        hlay_params.addLayout(vlay_labels1)
        hlay_params.addLayout(vlay_les1)
        hlay_params.addLayout(vlay_labels2)
        hlay_params.addLayout(vlay_les2)

        vlay.addLayout(hlay_params)
        
        hlay_btn = QtWidgets.QHBoxLayout()
        hlay_btn.addWidget(self.titcurve_btn)
        hlay_btn.addWidget(self.export_btn)
        vlay.addLayout(hlay_btn)

        vlay.addStretch(1)

        self.grid_container.addLayout(vlay, 0, 0)

    def connect_event_handlers(self):
        self.acid_combo.currentTextChanged.connect(self.read_param_le)
        self.acid_presets_combo.currentTextChanged.connect(self.update_k_values)
        for curr_name in self.params.index:
            if curr_name != 'acid_base':
                self.les[curr_name].textChanged.connect(self.read_param_le)

        self.titcurve_btn.clicked.connect(self.draw_titration_curve)

    def center(self):  # centers object on screen
        qr = self.frameGeometry()
        cp = QtWidgets.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def update_k_values(self):
        if self.acid_presets_combo.currentText() != '':
            for curr_k_lbl, curr_k in zip(self.params.index[1:4], k_values[self.acid_combo.currentText()][self.acid_presets_combo.currentText()]):
                self.les[curr_k_lbl].setText(str(curr_k))

    def read_param_le(self):
        for curr_name in self.params.index:
            if curr_name == 'acid_base':
                self.params[curr_name] = self.acid_combo.currentText()
            else:
                self.params[curr_name] = float(self.les[curr_name].text())

        self.titration.set_basic_params(
            self.params.iloc[1:4].tolist(), mode=self.params['acid_base'],
            kw=self.params.iloc[4])

    def draw_titration_curve(self):
        self.titration.curve(
            self.params['<i>V</i><sub>solution</sub> [L]'],
            self.params['<i>c</i><sub>solution</sub> [mol/L]'],
            self.params['<i>c</i><sub>titrant</sub> [mol/L]'],
            data_points=500, indep_var='pH',
            indep_var_min=self.params['pH<sub>min</sub>'],
            indep_var_max=self.params['pH<sub>max</sub>'])

        self.error_line.setData(
            self.titration.latest_curve[0],
            self.titration.latest_curve[1])

if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = titration_window(7.4)
    w.show()
    sys.exit(app.exec_())