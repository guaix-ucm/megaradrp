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

import math

import numpy as np
from scipy.stats import norm
import scipy.interpolate as ii
from scipy.ndimage.filters import convolve1d

from .vph import MegaraVPH
from .detector import MegaraDetectorSat
from .detector import ReadParams
from .focalplane import FocalPlane
from .fiberbundle import FiberBundle


class PseudoSlit(object):
    def __init__(self):

        # Positions of the fibers in the PS slit
        self.y_pos = {}

    def connect_fibers(self, fibid, pos):
        self.y_pos = dict(zip(fibid, pos))


class InternalOptics(object):
    def __init__(self):

        # Include the transmission of the fibers (all fibers are equal)
        #rawdata = np.loadtxt('tspect_0.1aa.dat')
        #self.trans_interp = ii.interp1d(rawdata[:,0] / 1e4, rawdata[:,1],
        #                                bounds_error=False, fill_value=0.0, copy=False)
        pass

    def transmission(self, wl):
        return np.ones_like(wl)
        #return self.trans_interp(wl)


class MegaraInstrument(object):
    def __init__(self):

        DSHAPE = (2056 * 2, 2048 * 2)
        PSCAN = 50
        OSCAN = 50

        layouttable = np.loadtxt('LCB_spaxel_centers.dat')

        fib_ids = layouttable[:,4].astype('int').tolist()
        bun_ids = layouttable[:,3].astype('int').tolist()

        self.fibers = FiberBundle(fib_ids, bun_ids)

        self.pseudo_slit = PseudoSlit()
        self.pseudo_slit.connect_fibers(fib_ids, layouttable[:,2])

        self.focal_plane = FocalPlane()
        self.focal_plane.connect_fibers(fib_ids, layouttable[:,0:2])

        self.vph = MegaraVPH()
        self.internal_optics = InternalOptics()

        self._internal_focus_factor = 1.0
        self._ref_focus = 123.123
        self._internal_focus = self._ref_focus

        eq = 1.0 * np.ones(DSHAPE)
        #eq[50:55,50:70] = 0.0

        dcurrent = 3.0 / 3600

        # eq = np.random.normal(loc=0.80, scale=0.01, size=(4096,4096))
        # eq = np.clip(eq, 0.0, 1.0)
        # fits.writeto('eq.fits', eq, clobber=True)
        # eq = fits.getdata('eq.fits')

        readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000)
        readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005)

        self.detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq, dark=dcurrent,
                                        readpars1=readpars1, readpars2=readpars2, bins='11')

    def set_cover(self, status):
        """Cover in the focal plane."""

        # This should be in other object
        self.focal_plane.set_cover(status)
        # Update visible fibers in the ps_slit

    def set_focus(self, x):
        """Arbitrary parametrization of the focus"""

        if  x < 118 or x > 128:
            raise ValueError('focus out of limits')

        self._internal_focus_factor = 1.0 + (math.cosh((x-self._ref_focus) / 2.0) - 1.0 ) / 5.0
        self._internal_focus = x

    def meta(self):
        _meta = {'detector': self.detector.metadata(),
                 'vph': self.vph.metadata(),
                 'focus': self._internal_focus,
                 'fplane': self.focal_plane.meta(),
                 }

        return _meta

    def project_rss(self, sigma, wl_in, spec_in):

        # This will compute only the illuminated fibers
        visible_fib_ids = [fibid for fibid, _ in self.focal_plane.get_visible_fibers()]
        return project_rss(visible_fib_ids, self.pseudo_slit, self.vph, self.detector, sigma, wl_in, spec_in)

    def project_rss_w(self):
        # This will compute only the illuminated fibers
        visible_fib_ids = [fibid for fibid, _ in self.focal_plane.get_visible_fibers()]
        sigma = self.fibers.sigma
        return project_rss_w(visible_fib_ids, self.pseudo_slit, self.vph, sigma)

    def simulate_focal_plane(self, wltable_in, photons_in):

        """Simulate the image in the focal plane,"""

        # Input spectra in each fiber, in photons
        # Fiber flux is 0 by default
        base_coverage = np.zeros((self.fibers.N, 1))

        for fibid, cover in self.focal_plane.get_visible_fibers():
            # FIBID starts in 1
            # Fill with covering factor
            # 0-0.5-1
            base_coverage[fibid-1] = cover

        allfiber_t = base_coverage * self.fibers.transmission(wltable_in)

        photons_in = allfiber_t * photons_in

        # Spectrograph optics transmission
        photons_in = photons_in * self.internal_optics.transmission(wltable_in)

        # VPH transmission
        photons_in = photons_in * self.vph.transmission(wltable_in)

        # QE of the detector
        photons_in = photons_in * self.detector.qe_wl(wltable_in)

        spec_in = photons_in

        return self.project_rss(self._internal_focus_factor* self.fibers.sigma, wltable_in, spec_in)


def project_rss(vis_fibs_id, pseudo_slit, vph, detector, sigma, wl_in, spec_in, scale=8):

    # Scale for super sampling
    # scale

    y_ps_fibers = [pseudo_slit.y_pos[fid] for fid in vis_fibs_id]
    DSHAPE = detector.dshape
    PIXSCALE = detector.pixscale
    # Scale for super sampling
    xcenter = detector.dshape[1] // 2
    ycenter = detector.dshape[0] // 2

    spos = -PIXSCALE * (np.arange(0, DSHAPE[1]*scale) - scale * xcenter) / scale

    wl_in_super = vph.ainv(y_ps_fibers, -spos, grid=True)
    wl_in_super = wl_in_super[:,::-1]

    spec_in_super = np.zeros_like(wl_in_super)

    # Resample to higher spatial resolution
    for i, fib_id in enumerate(vis_fibs_id):
        idx = fib_id - 1
        interpolator = ii.interp1d(wl_in, spec_in[idx])
        spec_in_super[i] = interpolator(wl_in_super[i])

    # kernel is constant in pixels
    # This is a gaussian convolved with a square
    kernel = compute_kernel(scale*sigma, truncate=5.0, d=0.5 * scale)
    out = convolve1d(spec_in_super, kernel, axis=1)

    # Downsample
    spec_in_detector = out[:,::scale]
    wl_in_detector = wl_in_super[:,::scale]

    ytrace = np.zeros_like(wl_in_detector)

    for idx, y_ps_fiber in enumerate(y_ps_fibers):
        ytrace[idx] = ycenter + vph.b(y_ps_fiber, wl_in_detector[idx]) / PIXSCALE

    nsig = 6 # At 6 sigma, the value of the profile * 60000 counts is
             # << 1
    final = np.zeros(DSHAPE)

    for idx, _ in enumerate(y_ps_fibers):
        minp = coor_to_pix(ytrace[idx].min() - nsig * sigma)
        maxp = coor_to_pix(ytrace[idx].max() + nsig * sigma)

        yp = np.arange(minp, maxp)

        base = pixcont_int_pix(yp[:,np.newaxis], ytrace[idx, :], sigma)

        #yspec = base * 7.5e3#spec_in_detector[idx]
        yspec = base * spec_in_detector[idx]
        final[minp:maxp,:] += yspec

    return final


def project_rss_w(visible_fib_ids, pseudo_slit, vph, detector, sigma):
    DSHAPE = detector.dshape
    PIXSCALE = detector.pixscale

    xcenter = DSHAPE[1] // 2
    ycenter = DSHAPE[0] // 2

    y_ps_fibers = [pseudo_slit.y_pos[fid] for fid in visible_fib_ids]

    spos = -PIXSCALE * (np.arange(0, DSHAPE[1]) - xcenter)
    wl_in_detector = vph.ainv(y_ps_fibers, -spos, grid=True)
    wl_in_detector = wl_in_detector[:,::-1]

    ytrace = np.zeros_like(wl_in_detector)

    for idx, y_pslit in enumerate(y_ps_fibers):
        ytrace[idx] = ycenter + vph.b(y_pslit, wl_in_detector[idx]) / PIXSCALE

    nsig = 6 # At 6 sigma, the value of the profile * 60000 counts is
             # << 1
    fvalue = 5e-6

    l1 = []
    l2 = []
    l3 = []
    l4 = []

    for idx, _ in enumerate(y_ps_fibers):
        minp = coor_to_pix(ytrace[idx].min() - nsig * sigma)
        maxp = coor_to_pix(ytrace[idx].max() + nsig * sigma)
        yp = np.arange(minp, maxp)
        base = pixcont_int_pix(yp[:,np.newaxis], ytrace[idx, :], sigma)
        sidx = np.nonzero(base >= fvalue)
        l1.extend(minp + sidx[0])
        l2.extend(sidx[1])
        l3.extend(base[sidx])
        l4.extend([idx]*len(sidx[0]))

    return l1, l2, l3, l4
            # row, col, val, fib

def coor_to_pix(x):
    return np.ceil(x -0.5).astype('int')


def pixcont_int_pix(i, x0, sig, d=0.5):
    """Integrate a gaussian profile."""
    zs = (i - x0) / sig
    z2 = zs + d / sig
    z1 = zs - d / sig
    return (norm.cdf(z2) - norm.cdf(z1)) / (2*d)


def pixcont_int(i, x0, sig):
    '''Integrate a gaussian profile.'''
    z2 = (i + 0.5 - x0) / sig
    z1 = (i - 0.5 - x0) / sig
    return norm.cdf(z2) - norm.cdf(z1)


def compute_kernel(sigma, d=0.5, truncate=5.0):
    """A kernel representing a Gaussian convolved with a square."""
    sd = float(sigma)
    # make the radius of the filter equal to truncate standard deviations
    lw = int(truncate * sd + 0.5)
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = pixcont_int(0, 0, sigma)
    sum = weights[lw]
    # calculate the kernel:
    for ii in range(1, lw + 1):
        tmp = pixcont_int_pix(ii, 0, sigma, d)
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum

    return weights
