import types

import astropy.units as u
from astropy.io import fits
import numpy
import pytest
from numina.exceptions import RecipeError

import megaradrp.recipes.scientific.base as base
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.recipes.calibration.lcbstdstar import LCBStandardRecipe
from megaradrp.testing.create_image import create_rss

NWL = 4300
EDGE = 150  # pixels without coverage at both ends
EXPTIME = 30.0
SKY = 50.0
STAR_MAG = 10.0
STAR_JY = (STAR_MAG * u.ABmag).to_value(u.Jy)
WL0, DWL = 6000.0, 0.1


def sensitivity(npix=NWL):
    """Counts / s / Jy, a smooth function of the pixel"""
    x = numpy.arange(npix)
    return 2.0e4 * (1.0 + 0.3 * (x - npix / 2) / npix)


@pytest.fixture(scope="module")
def lcb():
    hdrf = create_default_fiber_header("LCB")
    fp_conf = FocalPlaneConf.from_header(hdrf)
    center = min(fp_conf.connected_fibers(), key=lambda f: f.x**2 + f.y**2)
    return hdrf, fp_conf, center


def star_rss(lcb):
    """RSS with a standard star in the central fiber and a constant sky"""
    hdrf, fp_conf, center = lcb
    weights = numpy.zeros(623)
    for fiber in fp_conf.connected_fibers():
        # narrow enough to have all the flux in three rings
        dist2 = (fiber.x - center.x) ** 2 + (fiber.y - center.y) ** 2
        weights[fiber.fibid - 1] = numpy.exp(-dist2 / (2 * 0.3**2))
    weights /= weights.sum()
    star = EXPTIME * sensitivity() * STAR_JY
    data = weights[:, numpy.newaxis] * star + SKY
    wlmap = numpy.ones_like(data, dtype="int16")
    data[:, :EDGE] = 0
    data[:, -EDGE:] = 0
    wlmap[:, :EDGE] = 0
    wlmap[:, -EDGE:] = 0
    rss = create_rss(data.astype("float32"), wlmap)
    rss["FIBERS"].header = hdrf.copy()
    hdr = rss[0].header
    for key, val in [("CTYPE1", "AWAV"), ("CRPIX1", 1), ("CRVAL1", WL0), ("CDELT1", DWL), ("CUNIT1", "Angstrom")]:
        hdr[key] = val
    hdr["EXPTIME"] = EXPTIME
    hdr["AIRMASS"] = 1.0
    hdr["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    return rss


@pytest.fixture
def recipe(monkeypatch, lcb):
    reduced2d = fits.HDUList([fits.PrimaryHDU(numpy.zeros((10, 10)))])
    monkeypatch.setattr(base.ImageRecipe, "base_run", lambda self, rinput: (reduced2d, star_rss(lcb)))
    return LCBStandardRecipe()


def make_input(**kwds):
    values = dict(
        reference_spectrum=numpy.array([[5000.0, STAR_MAG], [8000.0, STAR_MAG]]),
        reference_spectrum_velocity=0.0,
        reference_extinction=numpy.array([[5000.0, 0.0], [8000.0, 0.0]]),
        sky_rss=None,
        ignored_sky_bundles=[],
        position=None,
        nrings=3,
        degrade_resolution_method="none",
        sigma_resolution=20.0,
        smoothing_knots=3,
    )
    values.update(kwds)
    return types.SimpleNamespace(**values)


def flux_valid(hdr):
    """0-based slice of the region of valid flux calibration"""
    return slice(hdr["PIXLIMF1"] - 1, hdr["PIXLIMF2"])


@pytest.mark.parametrize("method", ["none", pytest.param("fixed", marks=pytest.mark.slow)])
def test_lcb_standard(recipe, lcb, method):
    """The sensitivity of a star with known flux"""
    _, fp_conf, center = lcb
    result = recipe.run(make_input(degrade_resolution_method=method))

    # the 37 fibers around the star, the central one included
    assert len(result.fiber_ids) == 37
    assert center.fibid in result.fiber_ids
    assert result.sigma == (0 if method == "none" else 20.0)
    # the star without sky
    star = result.star_spectrum.open()[0].data
    expected_star = EXPTIME * sensitivity() * STAR_JY
    inner = slice(EDGE + 10, NWL - EDGE - 10)
    assert numpy.allclose(star[inner], expected_star[inner], rtol=2e-3)

    raw = result.sensitivity_raw.open()[0]
    region = flux_valid(raw.header)
    if method == "fixed":
        # the region of valid flux is at 2 sigma from the borders, where the
        # broadened spectrum is 2% lower; check at 4 sigma
        margin = int(2 * result.sigma)
        region = slice(region.start + margin, region.stop - margin)
    assert numpy.allclose(raw.data[region], sensitivity()[region], rtol=2e-3)
    master = result.master_sensitivity.open()[0]
    assert master.header["NUMTYPE"] == "MASTER_SENSITIVITY"
    assert master.header["TUNIT"] == "Jy"


@pytest.mark.slow
@pytest.mark.xfail(reason="the spline is fitted also where there is no coverage, with values near 0", strict=True)
def test_lcb_standard_smoothed_sensitivity(recipe):
    """The smoothed sensitivity follows the sensitivity in the region of valid flux"""
    result = recipe.run(make_input())
    master = result.master_sensitivity.open()[0]
    region = flux_valid(master.header)
    ratio = master.data[region] / sensitivity()[region]
    assert numpy.allclose(ratio, 1.0, rtol=0, atol=5e-3)


@pytest.mark.slow
def test_lcb_standard_position(recipe, lcb):
    """The star is extracted around the given position"""
    _, fp_conf, center = lcb
    result = recipe.run(make_input(position=(center.x, center.y)))
    assert center.fibid in result.fiber_ids


def test_lcb_standard_not_magnitudes(recipe):
    with pytest.raises(RecipeError):
        recipe.run(make_input(reference_spectrum=numpy.array([[5000.0, 1e3], [8000.0, 1e3]])))


def test_lcb_standard_bad_method(recipe):
    with pytest.raises(ValueError):
        recipe.run(make_input(degrade_resolution_method="other"))
