#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy
import scipy.interpolate as ii
import astropy.io.fits as fits
import astropy.wcs
from numina.instrument.simulation.efficiency import Efficiency


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


class EfficiencyFile(Efficiency):

    def __init__(self, fname, fill_value=0.0, factor=1.0):
        self.interpf = InterpolFile(fname, fill_value=fill_value, factor=factor)

    def response(self, wl):
        return self.interpf(wl)
