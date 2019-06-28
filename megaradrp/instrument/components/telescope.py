#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

from numina.instrument.hwdevice import HWDevice
from numina.instrument.simulation.efficiency import Efficiency

from megaradrp.instrument.components.instrument import FocusActuator


class M2FocusActuator(FocusActuator):
    def __init__(self, name, parent=None):
        super(M2FocusActuator, self).__init__(name, parent)

        # Focus
        self.internal_focus_factor = 1.0
        self._ref_focus = 1120
        self._internal_focus = self._ref_focus

    def set_focus(self, x):
        """Arbitrary parametrization of the focus"""
        if x < -4000 or x > 4000:
            raise ValueError('telescope focus out of limits')

        self.internal_focus_factor = 1+ 1.9 * (math.cosh((x - self._ref_focus) / 3000.0) - 1)
        self._internal_focus = x


class Telescope(HWDevice):

    def __init__(self, name, diameter, transmission=None):
        super(Telescope, self).__init__(name)

        self.diameter = diameter
        self.area = 0.25 * math.pi * diameter**2
        self.inc = None

        if transmission is None:
            self._transmission = Efficiency()
        else:
            self._transmission = transmission

        self.focus_actuator = M2FocusActuator('Focus', self)

    def transmission(self, wl):
        return self._transmission.response(wl)

    def connect(self, atmmodel):
        self.inc = atmmodel

    def set_focus(self, x):
        """Arbitrary parametrization of the focus"""

        self.focus_actuator.set_focus(x)