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
import scipy.interpolate as ii



class MegaraVPH(object):

    def __init__(self):
        self.SAMPLING = 9.0
        self.minwl = 3653.0
        self.maxwl = 4386.0
        self.wlmin_in = 0.3596
        self.wlmax_in = 0.4437
        self.res = 6028.0
        self.name = 'VPH405_LR'
        rr = np.loadtxt('VPH405_LR2.dat')
        r1 = rr[:,0] # Position in the pseudoslit
        r2 = rr[:,1] # WL
        r3 = rr[:,2] # X position
        r4 = rr[:,3] # Y position

        # Bivariate interpolations
        self.a = ii.SmoothBivariateSpline(r1, r2, r3)
        self.b = ii.SmoothBivariateSpline(r1, r2, r4)

        self.ainv = ii.SmoothBivariateSpline(r1, r3, r2)
        self.binv = ii.SmoothBivariateSpline(r1, r4, r2)

        # Include the transmission of the spectrograph
        #tvphraw = np.loadtxt('tvph_0.1aa.dat')
        #self.trans_interp = ii.interp1d(tvphraw[:,0] / 1e4, tvphraw[:,1],
        #                                bounds_error=False, fill_value=0.0,
        #                                copy=False)

    def distortion(self):
        pass

    def resolution(self, wl):
        # This is as VPH405_LR_res
        return  self.res * np.ones_like(wl)

    def metadata(self):
        return {'name': self.name}

    def wltable_interp(self):
        res_in = (self.wlmax_in / self.res) / self.SAMPLING
        return np.arange(self.wlmin_in, self.wlmax_in, res_in)

    def transmission(self, wl):
        return np.ones_like(wl)
        #return self.trans_interp(wl)
