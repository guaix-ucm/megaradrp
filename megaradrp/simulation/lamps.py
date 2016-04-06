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

import numpy as np
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda

from .wheel import HWDevice
from megaradrp.simulation.extended import create_th_ar_arc_spectrum


class Lamp(HWDevice):
    def __init__(self, name, factor=1.0, illumination=None):
        super(Lamp, self).__init__(name=name)
        self.factor = factor
        if illumination is None:
            self._illum = lambda x, y: np.ones_like(x)
        else:
            self._illum = illumination

    def flux(self, wl):
        return self.factor * np.ones_like(wl)

    def illumination(self, x, y):
        return self._illum(x, y)


class BlackBodyLamp(Lamp):

    def __init__(self, name, temp, factor=1.0, illumination=None):
        self.temp = temp
        super(BlackBodyLamp, self).__init__(name, factor=factor,
                                            illumination=illumination)

    def flux(self, wl_in):
        photons_in_flat = wl_in * blackbody_lambda(wl_in * u.um, self.temp) / 10.0
        return self.factor * photons_in_flat


class FlatLamp(Lamp):
    def __init__(self, name, photons, illumination=None):
        self.photons = photons
        super(FlatLamp, self).__init__(name, factor=photons,
                                       illumination=illumination)


class ArcLamp(Lamp):

    def flux(self, wl_in):
        return self.factor * create_th_ar_arc_spectrum(wl_in)
