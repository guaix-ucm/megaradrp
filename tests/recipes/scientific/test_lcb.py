import types

from astropy.io import fits
import numpy
import pytest
from numina.core import DataFrame
from numina.util.node import Node

import megaradrp.recipes.scientific.base as base
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.fluxcalib import PIXLIM_KEYS
from megaradrp.recipes.scientific.lcb import LCBImageRecipe
from megaradrp.testing.create_image import create_rss

NWL = 100
SKY = 5.0
OBJECT = 100.0
EXPTIME = 10.0
AIRMASS = 1.5


@pytest.fixture(scope="module")
def lcb():
    hdrf = create_default_fiber_header("LCB")
    return hdrf, FocalPlaneConf.from_header(hdrf)


@pytest.fixture(scope="module")
def sky_row(lcb):
    """Row of a sky fiber"""
    return lcb[1].sky_fibers()[0] - 1


def make_rss(lcb, sky=SKY):
    """RSS of the LCB, with a constant sky and an object in fiber 100"""
    hdrf, _ = lcb
    data = numpy.full((623, NWL), sky, dtype="float32")
    data[99] += OBJECT
    rss = create_rss(data, numpy.ones_like(data, dtype="int16"))
    rss["FIBERS"].header = hdrf.copy()
    hdr = rss[0].header
    for key, val in [("CTYPE1", "AWAV"), ("CRPIX1", 1), ("CRVAL1", 4000.0), ("CDELT1", 2.0), ("CUNIT1", "Angstrom")]:
        hdr[key] = val
    hdr["EXPTIME"] = EXPTIME
    hdr["AIRMASS"] = AIRMASS
    hdr["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    return rss


def make_sensitivity(value=4.0):
    hdu = fits.PrimaryHDU(numpy.full(NWL, value, dtype="float32"))
    hdu.header["UUID"] = "11111111-2222-3333-4444-555555555555"
    hdu.header["TUNIT"] = "Jy"
    for key, val in zip(PIXLIM_KEYS, [1, NWL, 1, NWL, 1, NWL]):
        hdu.header[key] = val
    return DataFrame(frame=fits.HDUList([hdu]))


@pytest.fixture
def recipe(monkeypatch, lcb):
    """The recipe, with the 2D reduction and the extraction replaced by a RSS"""
    reduced2d = fits.HDUList([fits.PrimaryHDU(numpy.zeros((10, 10)))])
    # the recipe calls the method of the base class with super()
    monkeypatch.setattr(base.ImageRecipe, "base_run", lambda self, rinput: (reduced2d, make_rss(lcb)))
    return LCBImageRecipe()


def make_input(**kwds):
    values = dict(ignored_sky_bundles=[], sky_rss=None, master_sensitivity=None, reference_extinction=None)
    values.update(kwds)
    return types.SimpleNamespace(**values)


def test_lcb_image(recipe, sky_row):
    """The sky of the SKY bundles is subtracted"""
    result = recipe.run(make_input())
    final = result.final_rss.open()[0].data
    reduced = result.reduced_rss.open()[0].data
    assert numpy.allclose(final[99], OBJECT)
    assert numpy.allclose(final[0], 0.0)
    assert numpy.allclose(reduced[99], OBJECT + SKY)
    # the sky RSS has the data of the sky fibers
    assert numpy.allclose(result.sky_rss.open()[0].data[sky_row], SKY)
    assert result.reduced_image.open()[0].data.shape == (10, 10)


def test_lcb_image_sky_rss(recipe, lcb):
    """The sky RSS is subtracted, if given"""
    sky_rss = make_rss(lcb, sky=2.0)
    sky_rss[0].data[99] = 2.0
    result = recipe.run(make_input(sky_rss=sky_rss))
    final = result.final_rss.open()[0].data
    assert numpy.allclose(final[99], OBJECT + SKY - 2.0)
    assert numpy.allclose(final[0], SKY - 2.0)


def test_lcb_image_flux_calibration(recipe, sky_row):
    """The three RSS are flux calibrated"""
    result = recipe.run(make_input(master_sensitivity=make_sensitivity(4.0)))
    for name, expected in [("final_rss", OBJECT), ("reduced_rss", OBJECT + SKY)]:
        img = getattr(result, name).open()
        assert numpy.allclose(img[0].data[99], expected / EXPTIME / 4.0)
        assert img[0].header["BUNIT"] == "Jy"
    assert numpy.allclose(result.sky_rss.open()[0].data[sky_row], SKY / EXPTIME / 4.0)


def test_lcb_image_extinction(recipe, sky_row):
    """final_rss and reduced_rss are corrected from extinction, not sky_rss"""
    # 0.1 mag / airmass at all wavelengths
    extinction = numpy.array([[3000.0, 0.1], [9000.0, 0.1]])
    result = recipe.run(make_input(reference_extinction=extinction))
    factor = 10 ** (0.4 * 0.1 * AIRMASS)
    assert numpy.allclose(result.final_rss.open()[0].data[99], OBJECT * factor)
    assert numpy.allclose(result.reduced_rss.open()[0].data[99], (OBJECT + SKY) * factor)
    assert numpy.allclose(result.sky_rss.open()[0].data[sky_row], SKY)


def test_base_run(monkeypatch):
    """The 2D reduction and the 1D reduction, with their intermediate results"""
    recipe = LCBImageRecipe()
    calls = []
    img2d = fits.HDUList([fits.PrimaryHDU(numpy.ones((10, 10)))])
    rss = fits.HDUList([fits.PrimaryHDU(numpy.ones((623, NWL)))])

    def combination(rinput, flow, method, method_kwargs):
        calls.append(("combination", method.__name__, method_kwargs))
        return img2d

    def reduction_1d(img, tracemap, wlcalib, fiberflat, twflat, offset):
        calls.append(("reduction_1d", tracemap, wlcalib, fiberflat, twflat, offset))
        return rss

    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: "flow")
    monkeypatch.setattr(base, "basic_processing_with_combination", combination)
    monkeypatch.setattr(recipe, "run_reduction_1d", reduction_1d)
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: calls.append(("save", name)))

    rinput = types.SimpleNamespace(
        obresult=types.SimpleNamespace(configuration=None),
        method="mean",
        method_kwargs={},
        master_apertures="traces",
        master_wlcalib="wlcalib",
        master_fiberflat="fiberflat",
        master_twilight=None,
        extraction_offset=[0.5],
    )
    reduced2d, reduced_rss = recipe.base_run(rinput)

    assert calls == [
        ("combination", "mean", {}),
        ("save", "reduced_image.fits"),
        ("reduction_1d", "traces", "wlcalib", "fiberflat", None, [0.5]),
        ("save", "reduced_rss.fits"),
    ]
    # the 2D image is a copy, with the headers of the recipe
    assert reduced2d is not img2d
    assert "NUMXVER" in reduced2d[0].header
    assert reduced_rss is rss


@pytest.mark.parametrize("twilight", [False, True])
def test_run_reduction_1d(monkeypatch, twilight):
    """The correctors of the 1D reduction, in order"""
    calls = []

    def corrector(name):
        class Corrector(Node):
            def __init__(self, *args, **kwds):
                super().__init__()

            def run(self, img):
                calls.append(name)
                return img

        return Corrector

    for name in ["ApertureExtractor", "FlipLR", "WavelengthCalibrator", "FiberFlatCorrector", "TwilightCorrector"]:
        monkeypatch.setattr(base, name, corrector(name))

    class Opener:
        def open(self):
            return "opened"

    recipe = LCBImageRecipe()
    img = fits.HDUList([fits.PrimaryHDU(numpy.ones((5, 5)))])
    result = recipe.run_reduction_1d(img, "traces", "wlcalib", Opener(), Opener() if twilight else None, offset=[0.0])
    expected = ["ApertureExtractor", "FlipLR", "WavelengthCalibrator", "FiberFlatCorrector"]
    if twilight:
        expected.append("TwilightCorrector")
    assert calls == expected
    assert result is img
