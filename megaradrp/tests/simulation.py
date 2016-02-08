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
import astropy.io.fits as fits


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
    def __init__(self, base, geom, directfun, readpars):

        self.base = base
        self.trim, self.pcol, self.ocol, self.orow = geom

        self.direcfun = directfun

        self.readpars = readpars

    def readout_in_buffer(self, elec, final):

        final[self.trim] = self.direcfun(elec[self.base])

        final[self.trim] = final[self.trim] / self.readpars.gain

        # We could use different RON and BIAS in each section
        for section in [self.trim, self.pcol, self.ocol, self.orow]:
            final[section] = self.readpars.bias + numpy.random.normal(final[section], self.readpars.ron)

        return final

class ReadParams(object):
    """Readout parameters of each channel."""
    def __init__(self, gain=1.0, ron=2.0, bias=1000.0):
        self.gain = gain
        self.ron = ron
        self.bias = bias


class MegaraDetector(object):
    """Simple MEGARA detector."""

    _binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
    _direc = ['normal', 'mirror']

    def __init__(self, shape, oscan, pscan, eq=1.0, dark=0.0, readpars1=None, readpars2=None,
                 bins='11', direction='normal'):

        if bins not in self._binning:
            raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

        if direction not in self._direc:
            raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

        if direction == 'normal':
            directfun = lambda x: x
        else:
            directfun = numpy.fliplr

        readpars1 = readpars1 if not None else ReadParams()
        readpars2 = readpars2 if not None else ReadParams()

        self.blocks = self._binning[bins]

        self.fshape, a0, geom1, geom2 = self.init_regions(shape, oscan, pscan, self.blocks)

        base1, base2 = a0

        self.virt1 = VirtualDetector(base1, geom1, directfun, readpars1)

        self.virt2 = VirtualDetector(base2, geom2, directfun, readpars2)

        self._det = numpy.zeros(shape, dtype='float64')

        self.eq = eq
        self.dark = dark
        # Exposure time since last reset
        self._time_last = 0.0

    def expose(self, source=0.0, time=0.0):
        self._time_last = time
        self._det += (source + self.dark) * time

    def readout(self):
        """Readout the detector."""

        elec_mean = self.eq * self._det
        elec = numpy.random.poisson(elec_mean)
        elec_pre = self.saturate(elec)
        # Do binning in the array
        elec_v = binning(elec_pre, self.blocks[0], self.blocks[1])
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

    def saturate(self, x):
        """Compute non-linearity and saturation."""
        a = -1.2 / 45000.0
        #b = 1.0
        #c = -7000.0
        #det = math.sqrt(b**2-4*a*c)
        #u1 = (-b+det) / (2*a)

        f1 = lambda x:x

        def f2(x):
            return x + a* (x - 45000)**2

        f3 = lambda x: 52000

        p1 = x < 45000
        p3 = x >= 54312 # 45000 + u1
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

        base1 = (slice(0, nr), slice(0,nc2))
        base2 = (slice(nr, nr2), slice(0,nc2))

        geom1 = trim1, pcol1, ocol1, orow1
        geom2 = trim2, pcol2, ocol2, orow2

        return (fshape, (base1, base2), geom1, geom2)

    def metadata(self):
        return {'exposed': self._time_last}


class MegaraImageFactory(object):
    CARDS_P = [
        ('OBSERVAT', 'ORM', 'Name of observatory'),
        ('TELESCOP', 'GTC', 'Telescope id.'),
        ('INSTRUME', 'MEGARA', 'Name of the Instrument'),
        ('ORIGIN', 'Simulator', 'FITS file originator'),
    ]

    def create(self, mode, meta, data):
        pheader = fits.Header(self.CARDS_P)

        # Detector
        meta_det = meta.get('detector', {})

        exptime = meta_det.get('exposed', 0.0)
        pheader['EXPTIME'] = exptime
        pheader['EXPOSED'] = exptime

        # VPH
        meta_vph = meta.get('vph', {})

        vph_name = meta_vph.get('name', 'unknown')
        pheader['VPH'] = vph_name

        # Focus
        focus = meta.get('focus', 0.0)
        pheader['FOCUS'] = focus

        # Focal plane
        meta_fplane = meta.get('fplane', {})
        cover = meta_fplane.get('cover', 'unknown')
        pheader['COVER'] = cover

        hdu1 = fits.PrimaryHDU(data, header=pheader)
        hdul = fits.HDUList([hdu1])
        return hdul


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


def simulate_bias_fits(factory, detector):
    """Simulate a FITS array."""
    final = simulate_bias(detector)

    meta = {'detector': detector.metadata()}
    data = final

    fitsfile = factory.create('bias', meta=meta, data=data)

    return fitsfile


def simulate_dark_fits(factory, detector, exposure):
    """Simulate a DARK FITS."""
    final = simulate_dark(detector, exposure)

    meta = {'detector': detector.metadata()}
    data = final

    fitsfile = factory.create('dark', meta=meta, data=data)

    return fitsfile


def simulate_flat_fits(factory, detector, exposure, source):
    """Simulate a FLAT FITS,"""
    final = simulate_flat(detector, exposure, source)

    meta = {'detector': detector.metadata()}
    data = final

    fitsfile = factory.create('flat', meta=meta, data=data)
    return fitsfile