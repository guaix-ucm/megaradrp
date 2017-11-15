#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

from astropy import units as u

from .efficiency import Efficiency


class LightFiber(object):
    def __init__(self, name, fibid, transmission=None, inactive=False):
        self.name = name
        self.fibid = fibid
        # Geometry of the fibers
        self.size = 0.31 * u.arcsec
        self.area = math.sqrt(3) * self.size ** 2 / 2.0
        self.fwhm = 3.6
        self.sigma = self.fwhm / 2.3548
        self.inactive = inactive

        if transmission is None:
            self._transmission = Efficiency()
        else:
            self._transmission = transmission

    def transmission(self, wl):
        return self._transmission.response(wl)

    def config_info(self):
        return {'name': self.name,
                'fibid': self.fibid,
                'inactive': self.inactive
                }


class FiberSet(object):
    def __init__(self, name, size, fwhm):
        self.fibers = {}
        self.bundles = {}
        self.name = name

        # Geometry of the fibers
        self.size = size
        self.area = math.sqrt(3) * self.size ** 2 / 2.0
        self.fwhm = fwhm
        self.sigma = self.fwhm / 2.3548

    @property
    def nfibers(self):
        return len(self.fibers)

    def transmission(self, wlin):
        # Loop over fibers
        import numpy

        result = numpy.zeros((self.nfibers, wlin.shape[0]))
        for lf in self.fibers.values():
            result[lf.fibid - 1] = lf.transmission(wlin)
        return result

    def config_info(self):
        return {'name': self.name,
                'nfibers': self.nfibers,
                'fwhm': self.fwhm,
                'sigma': self.sigma,
                }