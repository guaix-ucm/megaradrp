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

import numpy as np
from scipy.stats import norm
import scipy.interpolate as ii
from scipy.ndimage.filters import convolve1d

from .efficiency import Efficiency
from .device import HWDevice


class PseudoSlit(HWDevice):
    def __init__(self, name, insmode):

        super(PseudoSlit, self).__init__(name)

        # Positions of the fibers in the PS slit
        self.insmode = insmode
        self.fibers = {}

    def connect_fibers(self, fibers, positions):

        for fibid, pos in enumerate(positions, 1):
            lf = fibers[fibid]
            self.fibers[fibid] = (lf, pos)

    def y_pos(self, fibsid):
        result = []
        for fibid in fibsid:
            fiber, pos = self.fibers[fibid]
            result.append(pos)
        return result

    def config_info(self):
        return {'name': self.name, 'insmode': self.insmode}