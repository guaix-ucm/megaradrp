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
import astropy.io.fits as fits
import astropy.wcs


class InterpolFile(object):

    def __init__(self, fname, fill_value=0.0, factor=1.0):
        rawdata = numpy.loadtxt(fname)
        self._interp = ii.interp1d(rawdata[:,0] / 1e4, rawdata[:,1],
                                   bounds_error=False, fill_value=fill_value)
        self.factor = factor

    def __call__(self, wl):
        return self.factor * self._interp(wl)


class InterpolFitsUVES(object):
    """Interpolate spectrum in UVES format.

    This is the format of the sky spectrum file
    """
    def __init__(self, fname, fill_value=0.0):
        with fits.open(fname) as hdul:
            rawdata = hdul[0].data
            wcs = astropy.wcs.WCS(hdul[0].header)
            # WL is in Angstroms -> to microns
            all_wl = wcs.all_pix2world(range(hdul[0].data.size), 0)

        self._interp = ii.interp1d(all_wl[0] / 1e4, numpy.abs(rawdata),
                                   bounds_error=False,
                                   fill_value=fill_value)

    def __call__(self, wl):
        return self._interp(wl)


class Efficiency(object):

    def response(self, wl):
        return numpy.ones_like(wl)


class EfficiencyFile(Efficiency):

    def __init__(self, fname, fill_value=0.0, factor=1.0):
        self.interpf = InterpolFile(fname, fill_value=fill_value, factor=factor)

    def response(self, wl):
        return self.interpf(wl)
