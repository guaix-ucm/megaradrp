#
# Copyright 2016 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

import math

from .device import HWDevice
from .efficiency import Efficiency
from .instrument import FocusActuator


class M2FocusActuator(FocusActuator):
    def __init__(self, name, parent=None):
        super(M2FocusActuator, self).__init__(name, parent)

        # Focus
        self.internal_focus_factor = 1.0
        self._ref_focus = 0
        self._internal_focus = self._ref_focus

    def set_focus(self, x):
        """Arbitrary parametrization of the focus"""
        if x < -3000 or x > 3000:
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