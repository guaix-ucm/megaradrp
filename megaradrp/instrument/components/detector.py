#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

import numpy
from numpy.lib.stride_tricks import as_strided as ast
from numina.instrument.components.detector import DetectorBase, VirtualDetector


_logger = logging.getLogger(__name__)


def binning(arr, br, bc):
    """Return a binned view if 'arr'"""
    nr, nc = arr.shape

    if nr % br != 0 or nc % bc != 0:
        raise ValueError("'bin' must be an integer multiple of size")

    bnr = nr // br
    bnc = nc // bc
    m = arr.dtype.itemsize
    newshape = bnr, bnc, br, bc
    newstrides = nc*br*m, bc*m, nc*m, m
    binned = ast(arr, shape=newshape, strides=newstrides)
    return binned


class ReadParams(object):
    """Readout parameters of each channel."""
    def __init__(self, gain=1.0, ron=2.0, bias=1000.0):
        self.gain = gain
        self.ron = ron
        self.bias = bias


class MegaraDetector(DetectorBase):
    """Simple MEGARA detector."""

    _binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
    _direc = ['normal', 'mirror']

    def __init__(self, name, shape, oscan, pscan,
                 qe=1.0, qe_wl=None, dark=0.0,
                 readpars1=None, readpars2=None,
                 bins='11', direction='normal'):

        super(MegaraDetector, self).__init__(name, shape, qe, qe_wl, dark)

        self.oscan = oscan
        self.pscan = pscan

        self.readpars1 = readpars1 if not None else ReadParams()
        self.readpars2 = readpars2 if not None else ReadParams()

        self.bins = self._set_direction(direction)
        self._set_binning(bins)
        self.set_geometry()

    def configure(self, profile):
        if 'bins' in profile:
            self.set_binning(profile['bins'])

    def _set_binning(self, bins):

        if bins not in self._binning:
            raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

        self.bins = bins
        self.blocks = self._binning[self.bins]

    def set_binning(self, bins):

        self._set_binning(bins)
        self.set_geometry()

    def _set_direction(self, direction):

        if direction not in self._direc:
            raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

        if direction == 'normal':
            directfun = lambda x: x
        else:
            directfun = numpy.fliplr

        self.directfun = directfun

    def set_direction(self, direction):

        self._set_direction(direction)
        self.set_geometry()

    def set_geometry(self):
        self.fshape, (base1, base2), geom1, geom2 = self.init_regions(self.dshape,
                                                                      self.oscan, self.pscan,
                                                                      self.blocks)

        self.virt1 = VirtualDetector(base1, geom1, self.directfun, self.readpars1)
        self.virt2 = VirtualDetector(base2, geom2, self.directfun, self.readpars2)

    def pre_readout(self, elec_pre):
        # FIXME: there is a bug in numpy here
        # this crashes if you don't make a copy()
        elec_pre_inv = numpy.flipud(elec_pre).copy()
        elec_v = binning(elec_pre_inv, self.blocks[0], self.blocks[1])
        elec_p = elec_v.reshape(elec_v.shape[0], elec_v.shape[1], -1)
        elec_f = elec_p.sum(axis=-1)
        return elec_f

    def base_readout(self, elec_f):
        # Output image
        final = numpy.zeros(self.fshape)
        self.virt1.readout_in_buffer(elec_f, final)
        self.virt2.readout_in_buffer(elec_f, final)
        return final

    def saturate(self, x):
        """Compute non-linearity and saturation."""
        a = -1.2 / 45000.0
        # b = 1.0
        # c = -7000.0
        # det = math.sqrt(b**2-4*a*c)
        # u1 = (-b+det) / (2*a)

        f1 = lambda x:x

        def f2(xx):
            return xx + a* (xx - 45000)**2

        f3 = lambda x: 52000

        p1 = x < 45000
        p3 = x >= 54312  # 45000 + u1
        p2 = numpy.logical_and(~p1, ~p3)

        return numpy.piecewise(x, [p1, p2, p3], [f1, f2, f3])

    # @classmethod
    def init_regions(cls, detshape, oscan, pscan, bng):
        """Create a image with overscan for testing."""

        base_rows = detshape[0] // 2
        base_cols = detshape[1] // 2
        base_oscan = oscan
        base_pscan = pscan

        nr = base_rows // bng[0]
        nc = base_cols // bng[1]

        nr2 = 2 * nr
        nc2 = 2 * nc

        oscan1 = base_oscan // bng[0]
        oscan2 = oscan1 * 2

        psc1 = base_pscan // bng[0]
        psc2 = 2 * psc1
        fshape = (nr2 + oscan2, nc2 + psc2)

        # Row block 1
        rb1 = slice(0, nr)
        rb1m = slice(nr, nr + oscan1)

        # Row block 2
        rb2 = slice(nr + oscan2, nr2 + oscan2)
        rb2m = slice(nr + oscan1, nr + oscan2)

        # Col block
        cb = slice(psc1, nc2 + psc1)
        # Col block left
        cbl = slice(0, psc1)
        # Col block right
        cbr = slice(nc2 + psc1, nc2 + psc2)
        # Col block full
        cbf = slice(0, nc2 + psc2)

        # Mode normal
        trim1 = (rb1, cb)
        # prescan
        pcol1 = (rb1, cbl)
        # overscan
        ocol1 = (rb1, cbr)
        orow1 = (rb1m, cbf)

        trim2 = (rb2, cb)
        # prescan
        pcol2 = (rb2, cbr)
        # overscan
        ocol2 = (rb2, cbl)
        orow2 = (rb2m, cbf)

        base1 = (slice(0, nr), slice(0, nc2))
        base2 = (slice(nr, nr2), slice(0, nc2))

        geom1 = trim1, pcol1, ocol1, orow1
        geom2 = trim2, pcol2, ocol2, orow2

        return (fshape, (base1, base2), geom1, geom2)

    def init_config_info(self):
        info = super(MegaraDetector, self).init_config_info()
        info['exposed'] = self._time_last
        info['vbin'] = int(self.bins[-1])
        info['hbin'] = int(self.bins[0])
        return info


class MegaraDetectorSat(MegaraDetector):

    def saturate(self, x):
        y = super(MegaraDetectorSat, self).saturate(x)

        # Some pixels have a special nonlineary given by
        # y[100,3000] = self.esp_nonlinearity(x[100,3000])
        # y[3000,3000] = self.esp_nonlinearity(x[3000,3000])
        return y

    def esp_nonlinearity(self, x):
        sat = 12000.0
        return sat * (1-sat / (sat + x))

