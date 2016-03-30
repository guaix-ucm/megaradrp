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

import numpy
import scipy.interpolate as ii

from .efficiency import Efficiency

class MegaraVPH(object):

    def __init__(self, name, vphtable, resolution, transmission=None):
        self.SAMPLING = 9.0

        self.name = name
        self._res = resolution

        rr = numpy.loadtxt(vphtable)
        r1 = rr[:,0] # Position in the pseudoslit
        r2 = rr[:,1] # WL
        r3 = rr[:,2] # X position
        r4 = rr[:,3] # Y position

        self.wlmin = rr[:,1].min()
        self.wlmax = rr[:,1].max()

        # Bivariate interpolations
        self.ps_wl_x = ii.SmoothBivariateSpline(r1, r2, r3)
        self.ps_wl_y = ii.SmoothBivariateSpline(r1, r2, r4)

        self.ps_x_wl = ii.SmoothBivariateSpline(r1, r3, r2)
        self.ps_y_wl = ii.SmoothBivariateSpline(r1, r4, r2)

        # extreme values for interpolation

        self.wlmin_in = 0.98 * self.wlmin
        self.wlmax_in = 1.02 * self.wlmax

        if transmission is None:
            self._transmission = Efficiency()
        else:
            self._transmission = transmission

    def distortion(self):
        pass

    def resolution(self, wl):
        return  self._res.response(wl)

    def config_info(self):
        return {'name': self.name}

    def wltable_interp(self):
        res_in = (self.wlmax/ self.resolution(self.wlmax)) / self.SAMPLING
        return numpy.arange(self.wlmin_in, self.wlmax_in, res_in)

    def transmission(self, wl):
        return self._transmission.response(wl)


class DummyVPH(object):

    def __init__(self, name):
        self.SAMPLING = 9.0
        self.minwl = 3653.0
        self.maxwl = 4386.0
        self.wlmin_in = 0.3596
        self.wlmax_in = 0.4437
        self.res = 6028.0
        self.name = name


    def resolution(self, wl):
        # This is as VPH405_LR_res
        return  self.res * numpy.ones_like(wl)

    def config_info(self):
        return {'name': self.name}

    def wltable_interp(self):
        res_in = (self.wlmax_in / self.res) / self.SAMPLING
        return numpy.arange(self.wlmin_in, self.wlmax_in, res_in)

    def transmission(self, wl):
        return numpy.ones_like(wl)
