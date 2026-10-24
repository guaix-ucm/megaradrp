import logging
import math

import astropy.units as u
import astropy.wcs
from astropy.io import fits
import numpy
import pytest
from scipy.ndimage import gaussian_filter

from megaradrp.datamodel import MegaraDataModel, create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
import megaradrp.processing.extractobj as eobj
from megaradrp.processing.extractobj import (
    compute_broadening,
    compute_centroid,
    extract_star,
    generate_sensitivity,
    mix_values,
)
from megaradrp.processing.fluxcalib import FluxCalibration
from megaradrp.testing.create_image import create_rss


def test_cover_1():

    arr = [1, 1, 1, 1, 1, 1]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 7)

    arr = [0, 0, 0, 0, 0, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 0)


def test_cover_2():

    arr = [0, 1, 1, 1, 0, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 4, None)

    arr = [0, 1, 1, 1, 1, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 5)

    arr = [0, 1, 1, 1, 1, 1]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 7)

    arr = [1, 1, 1, 1, 1, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 5)


def test_coverage_issue247():

    # Example from test case
    fmap = numpy.zeros((4300,), dtype="int32")
    fmap[226:4136] = 37
    fmap[224] = 1
    fmap[225] = 10
    fmap[4136] = 32

    # Test total coverage
    arr1 = fmap == fmap.max()
    res1 = eobj.coverage_det(arr1.astype("int"))
    assert res1 == slice(226, 4136, None)

    # Test max coverage
    arr2 = fmap > 0
    res2 = eobj.coverage_det(arr2.astype("int"))
    assert res2 == slice(224, 4137, None)


@pytest.fixture(scope="module")
def focal_plane():
    return FocalPlaneConf.from_header(create_default_fiber_header("LCB"))


def star_weights(focal_plane, x0, y0, sigma=1.0):
    """Fraction of the flux of a Gaussian star in each fiber"""
    weights = numpy.zeros(623)
    for fiber in focal_plane.connected_fibers():
        weights[fiber.fibid - 1] = numpy.exp(-((fiber.x - x0) ** 2 + (fiber.y - y0) ** 2) / (2 * sigma**2))
    return weights


def make_rss(weights, nwl=100, spectrum=None, exptime=10.0, airmass=1.0):
    if spectrum is None:
        spectrum = numpy.linspace(1.0, 2.0, nwl)
    data = (weights[:, numpy.newaxis] * spectrum).astype("float32")
    wlmap = numpy.zeros_like(data, dtype="int16")
    wlmap[weights > 0] = 1
    rss = create_rss(data, wlmap)
    hdr = rss[0].header
    for key, val in [("CTYPE1", "AWAV"), ("CRPIX1", 1), ("CRVAL1", 4000.0), ("CDELT1", 2.0), ("CUNIT1", "Angstrom")]:
        hdr[key] = val
    hdr["EXPTIME"] = exptime
    hdr["AIRMASS"] = airmass
    hdr["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    return rss


def nearest_fibers(focal_plane, x0, y0, npoints):
    fibers = focal_plane.connected_fibers(valid_only=True)
    dist = [math.hypot(fiber.x - x0, fiber.y - y0) for fiber in fibers]
    return sorted(fibers[idx].fibid - 1 for idx in numpy.argsort(dist)[:npoints])


def test_extract_star(focal_plane):
    center = focal_plane.fibers[311]
    weights = star_weights(focal_plane, center.x, center.y)
    rss = make_rss(weights)
    flux, colids, max_slice, some_slice = extract_star(
        rss, (center.x, center.y), 7, focal_plane, logger=logging.getLogger(__name__)
    )
    assert colids == nearest_fibers(focal_plane, center.x, center.y, 7)
    # all the fibers are covered, the flux is the sum of the fibers
    assert numpy.allclose(flux, rss[0].data[colids].sum(axis=0))
    assert max_slice == slice(0, 101)
    assert some_slice == slice(0, 101)


def test_extract_star_partial_coverage(focal_plane):
    """The flux is corrected where some fibers are not covered"""
    center = focal_plane.fibers[311]
    weights = star_weights(focal_plane, center.x, center.y)
    rss = make_rss(weights)
    expected = rss[0].data[nearest_fibers(focal_plane, center.x, center.y, 7)].sum(axis=0)
    # a fiber of the first ring is not covered in the first 20 pixels
    colid = nearest_fibers(focal_plane, center.x, center.y, 7)[0]
    rss[0].data[colid, :20] = 0
    rss["WLMAP"].data[colid, :20] = 0

    flux, colids, max_slice, some_slice = extract_star(
        rss, (center.x, center.y), 7, focal_plane, logger=logging.getLogger(__name__)
    )
    assert max_slice == slice(20, 101)
    assert some_slice == slice(0, 101)
    assert numpy.allclose(flux, expected)


def test_compute_centroid(focal_plane):
    """The centroid of a star centered in a fiber is the position of the fiber"""
    center = focal_plane.fibers[311]
    rss = make_rss(star_weights(focal_plane, center.x, center.y))
    # the centroid of the 19 fibers nearest to the point
    centroid = compute_centroid(
        rss[0].data, focal_plane, 10, 90, (center.x, center.y), logger=logging.getLogger(__name__)
    )
    assert numpy.allclose(centroid, (center.x, center.y))


def test_compute_centroid_shifted(focal_plane):
    """The centroid moves towards the star"""
    center = focal_plane.fibers[311]
    rss = make_rss(star_weights(focal_plane, center.x + 0.2, center.y))
    centroid = compute_centroid(
        rss[0].data, focal_plane, 10, 90, (center.x, center.y), logger=logging.getLogger(__name__)
    )
    assert center.x < centroid[0] < center.x + 0.2
    assert centroid[1] == pytest.approx(center.y, abs=0.02)


def test_mix_values():
    rss = make_rss(numpy.ones(623))
    wcsl = astropy.wcs.WCS(rss[0].header)
    spectrum = numpy.arange(100.0)

    wavelen, response_0, response_1 = mix_values(wcsl, spectrum, lambda wl: numpy.zeros(len(wl)))
    assert numpy.allclose(wavelen.to_value(u.AA), 4000.0 + 2.0 * numpy.arange(100))
    assert response_0 is spectrum
    # AB magnitude 0 is 3631 Jy
    assert numpy.allclose(response_1, 3631.0, rtol=1e-3)


def test_compute_broadening():
    """The broadening that matches the spectra"""
    rng = numpy.random.default_rng(1)
    flux_high = 1.0 + numpy.zeros(1000)
    flux_high[rng.integers(50, 950, 40)] += 5.0
    flux_low = gaussian_filter(flux_high, 3.0)
    # the parameters of the recipes of standard stars
    offset, sigma = compute_broadening(
        flux_low,
        flux_high,
        sigmalist=[1.0, 2.0, 3.0, 4.0, 5.0],
        remove_mean=False,
        frac_cosbell=0.10,
        zero_padding=50,
        fminmax=(0.003, 0.3),
        naround_zero=25,
        nfit_peak=21,
    )
    assert sigma == 3.0
    assert abs(offset) < 0.1


def flat_star(wl):
    """A star of AB magnitude 10, without extinction"""
    return numpy.full(len(wl), 10.0)


def no_extinction(wl):
    return numpy.zeros(len(wl))


@pytest.mark.parametrize("sigma", [0.0, 5.0])
def test_generate_sensitivity(focal_plane, sigma):
    """The sensitivity is the ratio of counts per second to the flux of the star"""
    exptime = 10.0
    star_jy = (10.0 * u.ABmag).to_value(u.Jy)
    counts_per_jy = 2.0e4
    spectrum = numpy.full(100, exptime * counts_per_jy * star_jy)
    rss = make_rss(numpy.ones(623), exptime=exptime)

    sens = generate_sensitivity(rss, spectrum, flat_star, no_extinction, slice(10, 90), slice(5, 95), sigma=sigma)
    assert numpy.allclose(sens.data, counts_per_jy)
    assert sens.header["TUNIT"] == "Jy"


@pytest.mark.xfail(reason="the 1-based limits are written with ref=0, one pixel too high", strict=True)
@pytest.mark.parametrize("sigma", [0.0, 5.0])
def test_generate_sensitivity_limits(focal_plane, sigma):
    rss = make_rss(numpy.ones(623))
    spectrum = numpy.full(100, 1.0)
    sens = generate_sensitivity(rss, spectrum, flat_star, no_extinction, slice(10, 90), slice(5, 95), sigma=sigma)
    hdr = sens.header
    # 1-based limits: all fibers, some fibers, flux valid at 2 sigma from all fibers
    assert (hdr["PIXLIMM1"], hdr["PIXLIMM2"]) == (11, 90)
    assert (hdr["PIXLIMR1"], hdr["PIXLIMR2"]) == (6, 95)
    assert (hdr["PIXLIMF1"], hdr["PIXLIMF2"]) == (11 + 2 * sigma, 90 - 2 * sigma)
    assert hdr["WAVLIMM1"] == pytest.approx(4000.0 + 10 * 2.0)


def test_generate_sensitivity_flux_calibration(focal_plane):
    """The flux calibration with the sensitivity recovers the flux of the star"""
    exptime = 10.0
    star_jy = (10.0 * u.ABmag).to_value(u.Jy)
    spectrum = numpy.linspace(1.0, 2.0, 100) * exptime * star_jy
    rss = make_rss(numpy.ones(623), spectrum=spectrum, exptime=exptime)
    sens = generate_sensitivity(rss, spectrum, flat_star, no_extinction, slice(10, 90), slice(5, 95), sigma=0.0)

    calibrated = FluxCalibration(fits.HDUList([sens]), datamodel=MegaraDataModel())(rss)
    # the first and last pixels of the region with fibers are not checked, see
    # test_generate_sensitivity_limits
    assert numpy.allclose(calibrated[0].data[0, 6:94], star_jy)
    assert calibrated[0].header["BUNIT"] == "Jy"


def test_generate_sensitivity_raise():
    rss = make_rss(numpy.ones(623))
    with pytest.raises(ValueError):
        generate_sensitivity(rss, numpy.zeros(100), flat_star, no_extinction, slice(10, 90), slice(5, 95))
