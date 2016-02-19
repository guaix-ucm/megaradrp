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


class Efficiency(object):

    def response(self, wl):
        return np.ones_like(wl)


class EfficiencyFile(Efficiency):

    def __init__(self, fname, fill_value=0.0):
        rawdata = np.loadtxt(fname)
        self._interp = ii.interp1d(rawdata[:,0] / 1e4, rawdata[:,1],
                                   bounds_error=False, fill_value=fill_value)

    def response(self, wl):
        return self._interp(wl)
