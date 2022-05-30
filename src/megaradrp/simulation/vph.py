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
from astropy import units as u
from numina.instrument.simulation.efficiency import Efficiency


class MegaraVPH(object):

    def __init__(self, name, setup, vphtable, resolution, transmission=None, conf=None):
        self.SAMPLING = 9.0

        self.name = name
        self.setup = setup
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

        if conf is None:
            conf = {}

        self.wl_range = conf.get('wl_range', [0.0, 0.0, 0.0])

    def distortion(self):
        pass

    def resolution(self, wl):
        return  self._res.response(wl)

    def config_info(self):
        return {'name': self.name, 'setup': self.setup, 'wl_range': self.wl_range}

    def wltable_interp(self):
        res_in = (self.wlmax/ self.resolution(self.wlmax)) / self.SAMPLING
        return numpy.arange(self.wlmin_in, self.wlmax_in, res_in) * u.micron

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
