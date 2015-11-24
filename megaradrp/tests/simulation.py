#
# Copyright 2015 Universidad Complutense de Madrid
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

"""Simple monocromatic simulation"""


import numpy

from numpy.lib.stride_tricks import as_strided as ast


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


class VirtualDetector(object):
    """Each of the channels."""
    def __init__(self, base, geom, directfun, gain, bias, ron):

        self.base = base
        self.trim, self.pcol, self.ocol, self.orow = geom

        self.direcfun = directfun

        self.bias = bias
        self.gain = gain
        self.ron = ron

    def readout_in_buffer(self, elec, final):

        final[self.trim] = self.direcfun(elec[self.base])

        final[self.trim] = final[self.trim] / self.gain

        # We could use different RON and BIAS in each section
        for section in [self.trim, self.pcol, self.ocol, self.orow]:
            final[section] = self.bias + numpy.random.normal(final[section], self.ron)

        return final

class MegaraDetector(object):
    """Simple MEGARA detector."""

    _binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
    _direc = ['normal', 'mirror']

    def __init__(self, shape, oscan, pscan, eq=1.0, dark=0.0, gain=1.0, bias=100, ron=2.0,
                 bins='11', direction='normal'):

        if bins not in self._binning:
            raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

        if direction not in self._direc:
            raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

        if direction == 'normal':
            directfun = lambda x: x
        else:
            directfun = numpy.fliplr

        self.blocks = self._binning[bins]

        self.fshape, a0, geom1, geom2 = self.init_regions(shape, oscan, pscan, self.blocks)

        base1, base2 = a0

        self.virt1 = VirtualDetector(base1, geom1, directfun, gain, bias, ron)

        self.virt2 = VirtualDetector(base2, geom2, directfun, gain, bias, ron)

        self._det = numpy.zeros(shape, dtype='float64')

        self.eq = eq
        self.dark = dark

    def expose(self, source=0.0, time=0.0):
        self._det += (source + self.dark) * time

    def readout(self):
        """Readout the detector."""

        elec_mean = self.eq * self._det
        elec = numpy.random.poisson(elec_mean)
        # Do binning in the array
        elec_v = binning(elec, self.blocks[0], self.blocks[1])
        elec_p = elec_v.reshape(elec_v.shape[0], elec_v.shape[1], -1)
        elec_f = elec_p.sum(axis=-1)

        # Output image
        final = numpy.zeros(self.fshape)

        self.virt1.readout_in_buffer(elec_f, final)
        self.virt2.readout_in_buffer(elec_f, final)

        final = numpy.clip(final, 0, 2**16-1)
        self.reset()
        return final.astype('uint16')

    def reset(self):
        """Reset the detector."""
        self._det[:] = 0.0

    @classmethod
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

        oscan1 = base_oscan / bng[0]
        oscan2 = oscan1 * 2

        psc1 = base_pscan / bng[0]
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

        # Mode normal
        trim1 = (rb1, cb)
        pcol1 = (rb1, cbl)
        ocol1 = (rb1, cbr)
        orow1 = (rb1m, cb)

        trim2 = (rb2, cb)
        pcol2 = (rb2, cbr)
        ocol2 = (rb2, cbl)
        orow2 = (rb2m, cb)
        base1 = (slice(0, nr), slice(0,nc2))
        base2 = (slice(nr, nr2), slice(0,nc2))

        geom1 = trim1, pcol1, ocol1, orow1
        geom2 = trim2, pcol2, ocol2, orow2

        return (fshape, (base1, base2), geom1, geom2)


def simulate_bias(detector):
    """Simulate a BIAS array."""
    detector.expose(source=0.0, time=0.0)
    final = detector.readout()
    return final


def simulate_dark(detector, exposure):
    """Simulate a DARK array,"""
    detector.expose(source=0.0, time=exposure)
    final = detector.readout()
    return final


def simulate_flat(detector, exposure, source):
    """Simulate a FLAT array,"""
    detector.expose(source=source, time=exposure)
    final = detector.readout()
    return final
