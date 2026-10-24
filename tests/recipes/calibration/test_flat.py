import types

from astropy.io import fits
import numpy
import pytest
import numina.exceptions
from numina.util.node import Node

import megaradrp.recipes.calibration.flat as flat
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.recipes.calibration.flat import FiberFlatRecipe, _smoothing_window_check
from megaradrp.testing.create_image import create_rss

NWL = 4300


@pytest.fixture(scope="module")
def transmission():
    return numpy.random.default_rng(0).uniform(0.8, 1.2, 623)


def flat_rss(transmission, lamp, invalid=()):
    """Calibrated RSS of a flat: the spectrum of the lamp times the transmission of each fiber"""
    data = (transmission[:, numpy.newaxis] * lamp).astype("float32")
    rss = create_rss(data, numpy.ones_like(data, dtype="int16"))
    hdrf = create_default_fiber_header("LCB")
    for fibid in invalid:
        hdrf[f"FIB{fibid:03d}_V"] = False
    rss["FIBERS"].header = hdrf
    return rss


def expected_flat(transmission, invalid=()):
    valid = [idx for idx in range(623) if idx + 1 not in invalid]
    return valid, transmission / transmission[valid].mean()


def test_obtain_fiber_flat(transmission):
    """With a flat lamp, the fiber flat is the relative transmission of each fiber"""
    rss = flat_rss(transmission, numpy.ones(NWL))
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=1e-4)
    # the extensions are copied
    assert "WLMAP" in result
    assert result is not rss


def test_obtain_fiber_flat_no_coverage(transmission):
    """The pixels without coverage are 1"""
    rss = flat_rss(transmission, numpy.ones(NWL))
    rss["WLMAP"].data[5, :100] = 0
    rss[0].data[5, :100] = 0
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    assert numpy.all(result[0].data[5, :100] == 1.0)


@pytest.mark.xfail(reason="the average spectrum is divided by the coverage of all the fibers", strict=True)
def test_obtain_fiber_flat_invalid_fibers(transmission):
    """The invalid fibers do not change the normalization"""
    invalid = list(range(1, 31))
    rss = flat_rss(transmission, numpy.ones(NWL), invalid=invalid)
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission, invalid)
    assert numpy.allclose(result[0].data[valid], expected[valid, numpy.newaxis], rtol=1e-3)


@pytest.mark.xfail(reason="the spline with the default smoothing does not follow the spectrum of the lamp", strict=True)
def test_obtain_fiber_flat_lamp_shape(transmission):
    """The spectrum of the lamp is removed"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    rss = flat_rss(transmission, lamp)
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=1e-2)


@pytest.mark.parametrize("value, valid", [(31, True), (5, True), (3, False), (32, False)])
def test_smoothing_window_check(value, valid):
    if valid:
        assert _smoothing_window_check(value) == value
    else:
        with pytest.raises(numina.exceptions.ValidationError):
            _smoothing_window_check(value)


def test_process_flat2d(monkeypatch):
    recipe = FiberFlatRecipe()
    calls = []
    img = fits.HDUList([fits.PrimaryHDU(numpy.ones((10, 10)))])

    def combination(rinput, flow, method, method_kwargs):
        calls.append((flow, method.__name__, method_kwargs))
        return img

    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: "flow")
    monkeypatch.setattr(flat, "basic_processing_with_combination", combination)
    rinput = types.SimpleNamespace(
        obresult=types.SimpleNamespace(configuration=None), method="median", method_kwargs={"a": 1}
    )
    result = recipe.process_flat2d(rinput)
    assert calls == [("flow", "median", {"a": 1})]
    assert result is img
    assert "NUMXVER" in img[0].header


def test_fiber_flat_recipe(monkeypatch, transmission):
    """The recipe returns the 2D image and the RSS before the flip and the wavelength calibration"""
    calls = []
    img2d = fits.HDUList([fits.PrimaryHDU(numpy.full((10, 10), 7.0))])
    raw_rss = flat_rss(transmission, numpy.ones(NWL))
    raw_rss[0].data[:, :10] = 3.0

    class Extractor(Node):
        def __init__(self, tracemap, datamodel, offset):
            super().__init__()
            calls.append(("extractor", tracemap, offset))

        def run(self, img):
            return raw_rss

    class Calibrator(Node):
        def __init__(self, wlcalib, datamodel):
            super().__init__()
            calls.append(("calibrator", wlcalib))

        def run(self, img):
            calls.append(("calibrate", float(img[0].data[0, -1])))
            return img

    recipe = FiberFlatRecipe()
    saved = []
    monkeypatch.setattr(recipe, "process_flat2d", lambda rinput: img2d)
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: saved.append(name))
    monkeypatch.setattr(flat, "ApertureExtractor", Extractor)
    monkeypatch.setattr(flat, "WavelengthCalibrator", Calibrator)

    rinput = types.SimpleNamespace(
        master_apertures="traces", master_wlcalib="wlcalib", extraction_offset=[0.5], smoothing_window=31
    )
    result = recipe.run(rinput)

    assert saved == ["reduced_image.fits", "reduced_rss.fits"]
    # the wavelength calibration receives the flipped RSS
    assert calls == [("extractor", "traces", [0.5]), ("calibrator", "wlcalib"), ("calibrate", 3.0)]
    assert numpy.allclose(result.reduced_image.open()[0].data, 7.0)
    # the reduced RSS is the copy before the flip
    assert numpy.allclose(result.reduced_rss.open()[0].data[:, :10], 3.0)
    master = result.master_fiberflat.open()
    assert master[0].header["NUMTYPE"] == "MasterFiberFlat"
    assert master[0].data.shape == (623, NWL)
