import types

from astropy.io import fits
import numpy
import pytest
from numina.exceptions import ValidationError
from numina.util.node import Node

import megaradrp.recipes.calibration.twilight as twilight
from megaradrp.recipes.calibration.twilight import TwilightFiberFlatRecipe, pixel_2d_check, pixel_2d_check_or_none

NWL = 300


def test_pixel_2d_check():
    assert pixel_2d_check([1, 2]) == [1, 2]
    assert pixel_2d_check_or_none(None) is None
    assert pixel_2d_check_or_none([3, 4]) == [3, 4]
    with pytest.raises(ValidationError):
        pixel_2d_check([1, 2, 3])
    with pytest.raises(ValidationError):
        pixel_2d_check_or_none([1])


def test_combine_median_scaled():
    """The images are scaled to the level of the brightest one before the median (#341)"""
    base = numpy.random.default_rng(0).uniform(90, 110, size=(20, 30))
    arrays = [base * factor for factor in [0.25, 0.5, 1.0]]
    result = TwilightFiberFlatRecipe().combine_median_scaled(arrays, dtype="float64")
    assert numpy.allclose(result[0], base, rtol=1e-3)
    # number of images combined
    assert numpy.all(result[2] == 3)


def make_wlcalib(missing=(), error=()):
    return types.SimpleNamespace(missing_fibers=list(missing), error_fitting=list(error), total_fibers=623)


def test_good_ids_mask():
    mask = TwilightFiberFlatRecipe().good_ids_mask(make_wlcalib(missing=[5], error=[7, 8]))
    assert mask.shape == (623,)
    assert numpy.flatnonzero(~mask).tolist() == [4, 6, 7]


@pytest.mark.xfail(reason="good_ids_mask adds error_fitting to the missing_fibers of the calibration", strict=True)
def test_good_ids_mask_does_not_modify_the_calibration():
    calibration = make_wlcalib(missing=[5], error=[7, 8])
    TwilightFiberFlatRecipe().good_ids_mask(calibration)
    assert calibration.missing_fibers == [5]


def test_process_flat2d(monkeypatch):
    recipe = TwilightFiberFlatRecipe()
    calls = []
    img = fits.HDUList([fits.PrimaryHDU(numpy.ones((10, 10)))])

    def combination(rinput, flow, method):
        calls.append((flow, method))
        return img

    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: "flow")
    monkeypatch.setattr(twilight, "basic_processing_with_combination", combination)
    rinput = types.SimpleNamespace(obresult=types.SimpleNamespace(configuration=None))
    assert recipe.process_flat2d(rinput) is img
    assert calls == [("flow", recipe.combine_median_scaled)]
    assert "NUMXVER" in img[0].header


def test_run_reduction_1d(monkeypatch):
    calls = []

    def corrector(name):
        class Corrector(Node):
            def __init__(self, *args, **kwds):
                super().__init__()

            def run(self, img):
                calls.append(name)
                return img

        return Corrector

    for name in ["ApertureExtractor", "FlipLR", "WavelengthCalibrator", "FiberFlatCorrector"]:
        monkeypatch.setattr(twilight, name, corrector(name))

    fiberflat = types.SimpleNamespace(open=lambda: "opened")
    img = fits.HDUList([fits.PrimaryHDU(numpy.ones((5, 5)))])
    result = TwilightFiberFlatRecipe().run_reduction_1d(img, "traces", "wlcalib", fiberflat, offset=[0.0])
    assert calls == ["ApertureExtractor", "FlipLR", "WavelengthCalibrator", "FiberFlatCorrector"]
    assert result is img


@pytest.fixture
def twilight_rss():
    """RSS of a twilight: a continuum of level 2 plus a level per fiber"""
    level = numpy.linspace(50.0, 150.0, 623)
    data = numpy.full((623, NWL), 2.0) + level[:, numpy.newaxis]
    return level, fits.HDUList([fits.PrimaryHDU(data)])


@pytest.fixture
def recipe(monkeypatch, twilight_rss):
    recipe = TwilightFiberFlatRecipe()
    saved = []
    recipe.saved = saved
    img2d = fits.HDUList([fits.PrimaryHDU(numpy.full((10, 10), 7.0))])
    monkeypatch.setattr(recipe, "process_flat2d", lambda rinput: img2d)
    monkeypatch.setattr(recipe, "run_reduction_1d", lambda img, *args, **kwds: twilight_rss[1])
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: saved.append(name))
    return recipe


def make_input(continuum_region=None, missing=(5,)):
    return types.SimpleNamespace(
        master_apertures="traces",
        master_wlcalib=make_wlcalib(missing=missing),
        master_fiberflat=None,
        extraction_offset=[0.0],
        normalize_region=[100, 200],
        continuum_region=continuum_region,
    )


@pytest.mark.parametrize("continuum_region", [None, [100, 100]])
def test_twilight_recipe(recipe, twilight_rss, continuum_region):
    """Each fiber is its mean in the region, normalized to the mean of the good fibers"""
    level, rss = twilight_rss
    result = recipe.run(make_input(continuum_region=continuum_region))

    assert recipe.saved == ["reduced_image.fits", "reduced_rss.fits"]
    collapse = level + 2.0
    good = numpy.ones(623, dtype=bool)
    good[4] = False
    expected = collapse / collapse[good].mean()
    master = result.master_twilightflat.open()
    assert numpy.allclose(master[0].data, expected[:, numpy.newaxis])
    assert master[0].header["NUMTYPE"] == "MasterTwilightFlat"
    assert numpy.allclose(result.reduced_image.open()[0].data, 7.0)
    assert result.reduced_rss.open()[0].data.shape == (623, NWL)


def test_twilight_recipe_continuum(recipe, twilight_rss):
    """The mean of the continuum region is subtracted"""
    level, rss = twilight_rss
    rss[0].data[:, 250:280] = 2.0
    result = recipe.run(make_input(continuum_region=[250, 280]))
    good = numpy.ones(623, dtype=bool)
    good[4] = False
    expected = level / level[good].mean()
    assert numpy.allclose(result.master_twilightflat.open()[0].data, expected[:, numpy.newaxis])
