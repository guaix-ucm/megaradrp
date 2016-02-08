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

"""Extended multiwavelength simulation"""


from __future__ import print_function

import math

import numpy as np
import scipy.interpolate as ii
from scipy.stats import norm
from scipy.ndimage.filters import convolve1d

import megaradrp.tests.simulation as sim

# Real detector
DSHAPE = (2056 * 2, 2048 * 2)
#DSHAPE = (100 * 2, 100 * 2)
PSCAN = 50
OSCAN = 50

BINR = 1
BINC = 1

SAMPLING = 9.0
# These are microns -> pixels
PIXSCALE = 15.0e-3

SHAPE = DSHAPE[0] // BINR, DSHAPE[1] // BINC
HSHAPE = SHAPE[0] // 2, SHAPE[1] // 2


class MegaraDetectorSat(sim.MegaraDetector):

    def saturate(self, x):
        y = super(MegaraDetectorSat, self).saturate(x)

        # Some pixels have a special nonlineary given by
        #y[100,3000] = self.esp_nonlinearity(x[100,3000])
        #y[3000,3000] = self.esp_nonlinearity(x[3000,3000])
        return y

    def esp_nonlinearity(self, x):
        sat = 12000.0
        return sat * (1-sat / (sat + x))

    def qe_wl(self, wl):
        """QE per wavelenght."""
        return np.ones_like(wl)


def coor_to_pix(x):
    return np.ceil(x -0.5).astype('int')


def pixcont_int(i, x0, sig):
    '''Integrate a gaussian profile.'''
    z2 = (i + 0.5 - x0) / sig
    z1 = (i - 0.5 - x0) / sig
    return norm.cdf(z2) - norm.cdf(z1)


def pixcont_int_pix(i, x0, sig, d=0.5):
    """Integrate a gaussian profile."""
    zs = (i - x0) / sig
    z2 = zs + d / sig
    z1 = zs - d / sig
    return (norm.cdf(z2) - norm.cdf(z1)) / (2*d)


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


class MegaraVPH(object):

    def __init__(self):

        self.minwl = 3653.0
        self.maxwl = 4386.0
        self.wlmin_in = 0.3596
        self.wlmax_in = 0.4437
        self.res = 6028.0
        self.name = 'VPH405_LR'
        rr = np.loadtxt('VPH405_LR2.dat')
        r1 = rr[:,0] # Position in the pseudoslit
        r2 = rr[:,1] # WL
        r3 = rr[:,2] # X position
        r4 = rr[:,3] # Y position

        # Bivariate interpolations
        self.a = ii.SmoothBivariateSpline(r1, r2, r3)
        self.b = ii.SmoothBivariateSpline(r1, r2, r4)

        self.ainv = ii.SmoothBivariateSpline(r1, r3, r2)
        self.binv = ii.SmoothBivariateSpline(r1, r4, r2)

        # Include the transmission of the spectrograph
        #tvphraw = np.loadtxt('tvph_0.1aa.dat')
        #self.trans_interp = ii.interp1d(tvphraw[:,0] / 1e4, tvphraw[:,1],
        #                                bounds_error=False, fill_value=0.0,
        #                                copy=False)

    def distortion(self):
        pass

    def resolution(self, wl):
        # This is as VPH405_LR_res
        return  self.res * np.ones_like(wl)

    def metadata(self):
        return {'name': self.name}

    def wltable_interp(self):
        res_in = (self.wlmax_in / self.res) / SAMPLING
        return np.arange(self.wlmin_in, self.wlmax_in, res_in)

    def transmission(self, wl):
        return np.ones_like(wl)
        #return self.trans_interp(wl)


class PseudoSlit(object):
    def __init__(self):

        # Positions of the fibers in the PS slit
        self.y_pos = {}

    def connect_fibers(self, fibid, pos):
        self.y_pos = dict(zip(fibid, pos))


class FocalPlane(object):

    NAMES = {'UNSET': 3, 'LEFT': 2, 'RIGHT': 1, 'SET': 0}
    CODES = {3: 'UNSET', 2: 'LEFT', 1: 'RIGHT', 0: 'SET'}
    def __init__(self):

        self._cover_status = 3

        self._filter_s = lambda x: True
        self._cover_s = lambda x: 1.0

        self.fibid = []
        self.pos = []

    def set_cover(self, mode):
        """Cover in the focal plane."""

        if mode.upper() not in self.NAMES:
            raise ValueError('"%s" mode is not recognized')

        self._cover_status = self.NAMES[mode.upper()]

        self._cover_s = lambda x: 1.0

        if self._cover_status == 3:
            self._filter_s = lambda x: True
            self._cover_s = lambda x: 1.0
        elif self._cover_status == 2:
            self._filter_s = lambda pos: pos[0] >= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] > 0.0 else 0.5
        elif self._cover_status == 1:
            self._filter_s = lambda pos: pos[0] <= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] < 0.0 else 0.5
        else:
            self._filter_s = lambda pos: False
            self._cover_s = lambda pos: 0.0

    def connect_fibers(self, fibid, pos):
        self.fibid = fibid
        self.pos = pos

    def get_visible_fibers(self):
        p1 = [(fid, self._cover_s(pos)) for fid, pos in zip(self.fibid, self.pos) if self._filter_s(pos)]
        return p1

    def meta(self):
        return {'cover': self.CODES[self._cover_status]}


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


def project_rss_w(visible_fib_ids, pseudo_slit, vph, sigma):

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

def project_rss(vis_fibs_id, pseudo_slit, vph, detector, sigma, wl_in, spec_in, scale=8):

    # Scale for super sampling
    # scale

    y_ps_fibers = [pseudo_slit.y_pos[fid] for fid in vis_fibs_id]

    # Scale for super sampling
    xcenter = DSHAPE[1] // 2
    ycenter = DSHAPE[0] // 2

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


def project_rss_old(pseudo_slit, vph, sigma, wl_in, spec_in):
    """As it's done in the simulator."""
    # mm in the detector
    # X seems flipped L-R
    xcenter = DSHAPE[1] // 2
    ycenter = DSHAPE[0] // 2
    pos = -PIXSCALE * (np.arange(DSHAPE[1]) - xcenter)

    # Array with WL coordinates
    # pos must be increasing, so I revert
    wl_in_detector = vph.ainv(pseudo_slit.layouttable[:,2], -pos, grid=True)
    wl_in_detector = wl_in_detector[:,::-1]

    spec_in_detector = np.zeros_like(wl_in_detector)
    for i in range(pseudo_slit.layouttable.shape[0]):
        interpolator = ii.interp1d(wl_in, spec_in[i])
        spec_in_detector[i] = interpolator(wl_in_detector[i])

    ytrace = np.zeros_like(wl_in_detector)
    for idx, fiber in enumerate(pseudo_slit.layouttable):
        ytrace[idx] = ycenter + vph.b(fiber[2], wl_in_detector[idx]) / PIXSCALE

    all = np.arange(DSHAPE[1])
    implimg = np.zeros(DSHAPE)

    # map the trace in two pixels
    pix1s = coor_to_pix(ytrace)
    fracs = ytrace - pix1s + 0.5
    implimg[pix1s, all] = fracs * spec_in_detector
    implimg[pix1s - 1, all] = (1 - fracs) * spec_in_detector

    # This is a gaussian convolved with a square
    kernel = compute_kernel(sigma, truncate=5.0)

    out = convolve1d(implimg, kernel, axis=0)

    return out


def create_dummy_spectrum(wl_in):
    wlc = [0.37, 0.38, 0.39, 0.40, 0.41, 0.42]
    s = 8e-6

    photons_in = np.zeros_like(wl_in)

    for c in wlc:
        photons_in +=  4e4 * np.exp(-0.5*((wl_in-c) / s)**2)

    photons_in += 4.0e4
    return photons_in


# A ThAR arc...
def create_th_ar_arc_spectrum(wl_in):

    lines = [(3719.41400, 58.52137),
             (3803.10800, 54.91138),
             (3828.37100, 56.63661),
             (3839.72400, 46.62416),
             (4019.13100, 40.06343),
             (4071.99600, 42.12776),
             (4131.76200, 50.57991),
             (4158.58400, 68.03227),
             (4200.65300, 48.00637),
             (4277.55800, 78.81951),
             (4348.11900, 80.24873)
             ]

    s = 8e-6
    #wl_in = vph.wltable_interp()
    photons_in = 0.0 * wl_in

    for cwl, flux in lines:
        c = cwl / 1e4
        photons_in +=  flux * 1.5e4 * np.exp(-0.5*((wl_in-c) / s)**2)

    #photons_in += 1.0e2

    return photons_in


class MegaraInstrument(object):
    def __init__(self):

        layouttable = np.loadtxt('LCB_spaxel_centers.dat')

        fib_ids = layouttable[:,4].astype('int').tolist()
        bun_ids = layouttable[:,3].astype('int').tolist()

        self.fibers = FiberBundle(fib_ids, bun_ids)

        self.pseudo_slit = PseudoSlit()
        self.pseudo_slit.connect_fibers(fib_ids, layouttable[:,2])

        self.focal_plane = FocalPlane()
        self.focal_plane.connect_fibers(fib_ids, layouttable[:,0:2])

        self.factory = sim.MegaraImageFactory()
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

        readpars1 = sim.ReadParams(gain=1.0, ron=2.0, bias=1000)
        readpars2 = sim.ReadParams(gain=1.0, ron=2.0, bias=1005)

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

    def simulate_flat(self, wltable, photons_in, exposure, repeat=1):
        """Simulate a FLAT/ARC"""

        out = self.simulate_focal_plane(wltable, photons_in)

        for i in range(repeat):
            self.detector.expose(source=out, time=exposure)

            final = self.detector.readout()

            fitsfile = self.factory.create('fiber-flat', meta=self.meta(), data=final)
            yield fitsfile

    def simulate_arc(self, wltable, photons_in, exposure, repeat=1):
        """Simulate a FLAT/ARC"""

        out = self.simulate_focal_plane(wltable, photons_in)

        for i in range(repeat):
            self.detector.expose(source=out, time=exposure)

            final = self.detector.readout()

            fitsfile = self.factory.create('fiber-flat', meta=self.meta(), data=final)
            yield fitsfile

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
