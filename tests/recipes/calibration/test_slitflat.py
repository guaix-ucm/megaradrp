import types

from astropy.io import fits
import numpy
import pytest

import megaradrp.recipes.calibration.slitflat as slitflat
from megaradrp.recipes.calibration.slitflat import SlitFlatRecipe

SHAPE = (200, 150)


@pytest.fixture
def pixel_qe():
    """Pixel to pixel variations of the efficiency, 2%"""
    return 1.0 + 0.02 * numpy.random.default_rng(3).standard_normal(SHAPE)


def run_recipe(monkeypatch, data):
    recipe = SlitFlatRecipe()
    calls = []
    img = fits.HDUList([fits.PrimaryHDU(data.astype("float32"))])
    img[0].header["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"

    def combination(rinput, flow, method):
        calls.append(("combination", flow, method.__name__))
        return img

    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: "flow")
    monkeypatch.setattr(slitflat, "basic_processing_with_combination", combination)
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: calls.append(("save", name)))
    monkeypatch.setattr(recipe, "save_intermediate_array", lambda arr, name: calls.append(("save", name)))
    rinput = types.SimpleNamespace(obresult=types.SimpleNamespace(configuration=None))
    return recipe.run(rinput), calls


def test_slitflat_recipe(monkeypatch, pixel_qe):
    """The master slit flat contains the pixel to pixel variations, not the illumination"""
    yy, xx = numpy.mgrid[: SHAPE[0], : SHAPE[1]]
    illumination = 1000.0 * (1 + 0.3 * xx / SHAPE[1]) * (1 + 0.2 * numpy.sin(yy / 40.0))
    result, calls = run_recipe(monkeypatch, illumination * pixel_qe)

    assert calls[0] == ("combination", "flow", "median")
    assert [name for kind, name in calls[1:]] == [
        "reduced_image.fits",
        "median_image1.fits",
        "qe_filter1.fits",
        "median_image2.fits",
        "qe_filter2.fits",
    ]
    master = result.master_slitflat.open()[0]
    assert master.header["NUMTYPE"] == "MasterSlitFlat"
    assert master.header["UUID"] != "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    assert master.data.dtype == numpy.float32
    # far from the borders of the image and of the channels
    inner = (slice(20, 80), slice(20, 130))
    ratio = master.data[inner] / pixel_qe[inner]
    assert numpy.median(ratio) == pytest.approx(1.0, abs=2e-3)
    assert numpy.std(ratio) < 0.01
    assert numpy.allclose(result.reduced_image.open()[0].data, illumination * pixel_qe, rtol=1e-5)


def test_slitflat_recipe_zeros(monkeypatch, pixel_qe):
    """The pixels without signal are 1, without warnings (#360)"""
    data = 1000.0 * pixel_qe
    data[:, :40] = 0.0
    result, _ = run_recipe(monkeypatch, data)
    assert numpy.all(result.master_slitflat.open()[0].data[:, :30] == 1.0)
