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

class FiberBundle(object):
    def __init__(self, fid, bid):
        # Geometry of the fibers
        self.size = 0.31
        self.fwhm = 3.6
        self.sigma = self.fwhm / 2.3548

        self.fibs_id = fid
        self.bund_id = bid

        self.N = len(self.fibs_id)
        # Include the transmission of the fibers (all fibers are equal)
        #tfiberraw = np.loadtxt('tfiber_0.1aa_20m.dat')
        #self.tfiber_interp = ii.interp1d(tfiberraw[:,0] / 1e4, tfiberraw[:,1], bounds_error=False, fill_value=0.0, copy=False)

    def transmission(self, wl):
        return np.ones_like(wl)
        #return self.tfiber_interp(wl)
