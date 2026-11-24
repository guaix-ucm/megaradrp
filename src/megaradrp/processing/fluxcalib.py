#
# Copyright 2017-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import datetime
import logging
import sys

if sys.version_info[:2] <= (3, 10):
    datetime.UTC = datetime.timezone.utc

import astropy.wcs
import astropy.units as u

import numpy
from numina.processing import Corrector

_logger = logging.getLogger(__name__)


PIXLIM_KEYS = ["PIXLIMF1", "PIXLIMF2", "PIXLIMR1", "PIXLIMR2", "PIXLIMM1", "PIXLIMM2"]
WAVLIM_KEYS = ["WAVLIMF1", "WAVLIMF2", "WAVLIMR1", "WAVLIMR2", "WAVLIMM1", "WAVLIMM2"]


class FluxCalibration(Corrector):
    """A Node that calibrates absolute flux."""

    def __init__(self, sensitivity, datamodel=None, dtype="float32"):

        self.sensp = sensitivity["primary"]
        self.pixlims = {key: (self.sensp.header[key] - 1) for key in PIXLIM_KEYS}

        calibid = datamodel.get_imgid(sensitivity)
        self.target_unit = self.sensp.header["TUNIT"]
        super(FluxCalibration, self).__init__(datamodel=datamodel, calibid=calibid, dtype=dtype)

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug("calibrate flux in image %s", imgid)
        hdr = img["primary"].header
        exptime = hdr["EXPTIME"]
        img["primary"].data /= exptime
        newdata = numpy.zeros_like(img["primary"].data)
        limr1 = self.pixlims["PIXLIMR1"]
        limr2 = self.pixlims["PIXLIMR2"]
        validr = slice(limr1, limr2 + 1, 1)
        # with numpy.errstate(invalid='ignore', divide='ignore'):
        newdata[:, validr] = img["primary"].data[:, validr] / self.sensp.data[validr]
        img["primary"].data = newdata
        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr["BUNIT"] = self.target_unit
        hdr["NUM-STD"] = self.calibid
        hdr["NUM-FLUX"] = self.calibid

        update_flux_limits(hdr, self.pixlims, ref=0)

        hdr["history"] = f"Flux calibration in {imgid}"
        tnow = datetime.datetime.now(datetime.UTC)
        hdr["history"] = f"Flux calibration time {tnow.isoformat()}"


def update_flux_limits(header, pixlims, wcs=None, ref=1):
    """Update keywords used for flux limits"""

    pixlim_coms = [
        "Start of valid flux calibration",
        "End of valid flux calibration",
        "Start of region with at least one fiber",
        "End of region with at least one fiber",
        "Start of region with all fibers",
        "End of region with all fibers",
    ]

    if ref not in [0, 1]:
        raise ValueError("ref must be 0 or 1")

    off = (ref + 1) % 2

    if wcs is None:
        wcs = astropy.wcs.WCS(header)

    for key, com in zip(PIXLIM_KEYS, pixlim_coms):
        header[key] = (pixlims[key] + off, com)

    r1 = numpy.array([pixlims[key] for key in PIXLIM_KEYS])
    r2 = numpy.zeros_like(r1)
    lm = numpy.array([r1, r2])
    # Values are 0-based
    wavelen_ = wcs.all_pix2world(lm.T, ref)
    if wcs.wcs.cunit[0] == u.dimensionless_unscaled:
        # CUNIT is empty, assume Angstroms
        wavelen = wavelen_[:, 0] * u.AA
    else:
        wavelen = wavelen_[:, 0] * wcs.wcs.cunit[0]

    for idx, (key, com) in enumerate(zip(WAVLIM_KEYS, pixlim_coms)):
        header[key] = (wavelen[idx].to(u.AA).value, com)

    return header


def windows_check(value):
    """Check that `value` is a list of wavelength windows [[wl1, wl2], ...], with wl1 < wl2"""
    from numina.exceptions import ValidationError

    for window in value:
        try:
            w1, w2 = (float(v) for v in window)
        except (TypeError, ValueError):
            raise ValidationError(f"each window must be a pair of numbers [wl1, wl2], got {window!r}")
        if not w1 < w2:
            raise ValidationError(f"the window {window!r} must have wl1 < wl2")
    return value


def sensitivity_skip_windows():
    """Default windows (air, Angstrom) skipped by the smoothing of the sensitivity

    The telluric bands of the atmosphere and the strong lines of the standard
    stars, from ``instrument/configs/sensitivity_skip_windows.ecsv``.
    """
    import pkgutil

    from astropy.table import Table

    data = pkgutil.get_data("megaradrp.instrument.configs", "sensitivity_skip_windows.ecsv")
    table = Table.read(data.decode("utf8"), format="ascii.ecsv")
    return [(float(a), float(b)) for a, b in zip(table["wl_min"], table["wl_max"])]


def _supported(spline, x, weights, min_points=2):
    """The fit is finite and each interval between knots has pixels with weight"""
    if not numpy.all(numpy.isfinite(spline(x))):
        return False
    edges = spline.get_knots()
    used = x[weights > 0]
    counts, _ = numpy.histogram(used, bins=edges)
    return bool(numpy.all(counts >= min_points))


def smooth_sensitivity(wl, sens_raw, header, knots=5, adaptive=False, skip_windows=(), min_fraction=0.3):
    """The sensitivity smoothed with a spline fitted in the region of valid flux

    The spline is fitted only between ``PIXLIMF1`` and ``PIXLIMF2`` of
    `header`: outside it the raw sensitivity is 0 or affected by the partial
    coverage of the fibers, and fitting it biases the result (#357). Inside,
    the pixels in `skip_windows` (telluric bands, lines of the star) and those
    without a positive finite value have weight 0. Outside the region, the
    smoothed sensitivity continues with straight lines with the value and
    the slope of the spline at each border.

    If the windows leave less than `min_fraction` of the pixels of the
    region, they are ignored, with a warning.

    Parameters
    ----------
    wl : numpy.ndarray
        Wavelength of each pixel, in Angstrom
    sens_raw : numpy.ndarray
        Raw sensitivity
    header : astropy.io.fits.Header
        Header with the keywords PIXLIMF1 and PIXLIMF2 (1-based)
    knots : int or list of float
        Number of interior knots, or their wavelengths
    adaptive : bool
        Move the knots to their best positions; if False they are fixed
        (equidistant in the region if `knots` is a number)
    skip_windows : iterable of (float, float)
        Wavelength windows with weight 0

    Returns
    -------
    smoothed : numpy.ndarray
    knots : numpy.ndarray
        The interior knots used
    """
    from numina.array.numsplines import AdaptiveLSQUnivariateSpline

    wl = numpy.asarray(wl, dtype=float)
    sens_raw = numpy.asarray(sens_raw, dtype=float)
    start, stop = header["PIXLIMF1"] - 1, header["PIXLIMF2"]
    x, y = wl[start:stop], sens_raw[start:stop]
    usable = numpy.isfinite(y) & (y > 0)
    weights = usable.astype(float)
    skip = numpy.zeros_like(usable)
    for w1, w2 in skip_windows or ():
        skip |= (x >= w1) & (x <= w2)
    if numpy.count_nonzero(usable & ~skip) >= min_fraction * x.size:
        weights[skip] = 0.0
    else:
        _logger.warning(
            "the windows skipped by the smoothing leave %d of %d pixels, they are ignored",
            numpy.count_nonzero(usable & ~skip),
            x.size,
        )
    if not numpy.isscalar(knots):
        inside = [k for k in knots if x[0] < k < x[-1]]
        if len(inside) < len(knots):
            _logger.warning("knots outside the region of valid flux %.1f-%.1f are not used", x[0], x[-1])
        knots = inside
    ysafe = numpy.where(usable, y, 0.0)
    spline = AdaptiveLSQUnivariateSpline(x=x, y=ysafe, t=knots, w=weights, adaptive=adaptive)
    if numpy.any(skip & usable) and not _supported(spline, x, weights):
        # an interval between knots without data, inside a window
        _logger.warning("the windows leave knots of the smoothing without data, repeating without them")
        weights = usable.astype(float)
        spline = AdaptiveLSQUnivariateSpline(x=x, y=ysafe, t=knots, w=weights, adaptive=adaptive)
    smoothed = spline(wl)
    deriv = spline.derivative()
    for border, side in ((x[0], wl < x[0]), (x[-1], wl > x[-1])):
        if numpy.any(side):
            smoothed[side] = spline(border) + deriv(border) * (wl[side] - border)
    return smoothed, spline.get_knots()[1:-1]
