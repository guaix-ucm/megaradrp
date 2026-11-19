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


def test_good_ids_mask_contents():
    """A fiber without a solution is bad, even if the calibration does not list it as missing"""
    calibration = make_wlcalib(missing=[5])
    calibration.contents = [types.SimpleNamespace(fibid=fibid) for fibid in range(1, 623) if fibid != 5]
    mask = TwilightFiberFlatRecipe().good_ids_mask(calibration)
    assert numpy.flatnonzero(~mask).tolist() == [4, 622]


def test_good_ids_mask_does_not_modify_the_calibration():
    """The calibration is not modified (#360)"""
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


def test_select_normalize_region():
    """The initial window, moved out of the zero order in the direction in which the signal grows"""
    from megaradrp.recipes.calibration.twilight import select_normalize_region

    npix = 1000
    x = numpy.arange(npix)
    # the signal grows to the right
    data = numpy.tile(1.0 + x / npix, (10, 1))
    coverage = numpy.ones(data.shape, dtype=bool)
    zero = numpy.zeros(data.shape, dtype=bool)
    good = numpy.ones(10, dtype=bool)
    region = (400, 500)
    # without zero order, the initial window
    assert select_normalize_region(data, coverage, zero, good, region=region) == region
    # zero order in one good fiber: the window moves to the right, just after it
    zero[2, 450:470] = True
    assert select_normalize_region(data, coverage, zero, good, region=region) == (470, 570)
    # the signal grows to the left: the window moves to the left, just before it
    assert select_normalize_region(data[:, ::-1], coverage, zero, good, region=region) == (350, 450)
    # zero order in a bad fiber does not move the window
    good[2] = False
    assert select_normalize_region(data, coverage, zero, good, region=region) == region
    good[2] = True
    # no free window to the right: to the left
    zero[2, 450:470] = False
    zero[2, 450:] = True
    assert select_normalize_region(data, coverage, zero, good, region=region) == (350, 450)
    # a good fiber without coverage is like the zero order
    zero[:] = False
    coverage[3, 480:520] = False
    assert select_normalize_region(data, coverage, zero, good, region=region) == (520, 620)
    # no window at all
    coverage[3] = False
    assert select_normalize_region(data, coverage, zero, good, region=region) is None
    assert select_normalize_region(data, coverage, zero, good, region=(0, 2000)) is None


def test_masked_mean():
    """The masked pixels are left out without bias from the shape of the spectrum"""
    from megaradrp.recipes.calibration.twilight import masked_mean

    spectrum = numpy.array([1.0, 2.0, 4.0, 2.0])
    data = numpy.array([1.0, 2.0, 3.0])[:, numpy.newaxis] * spectrum
    good = numpy.ones(3, dtype=bool)
    usable = numpy.ones(data.shape, dtype=bool)
    # without masked pixels, the mean of each fiber
    assert numpy.allclose(masked_mean(data, usable, (0, 4), good), data.mean(axis=1))
    # the brightest pixel of the second fiber masked, and extra light there
    data[1, 2] += 10.0
    usable[1, 2] = False
    assert numpy.allclose(masked_mean(data, usable, (0, 4), good), [2.25, 4.5, 6.75])
    # a fiber without usable pixels gets its mean
    usable[2] = False
    assert numpy.isclose(masked_mean(data, usable, (0, 4), good)[2], data[2].mean())


def test_twilight_recipe_masked_automatic(recipe, twilight_rss):
    """With the automatic region, the window leaves the zero order towards more signal; invalid fibers left out"""
    from megaradrp.processing.fiberflat import flat_mask_hdu

    level, rss = twilight_rss
    x = numpy.arange(NWL)
    shape = 1.0 + numpy.exp(-0.5 * ((x - 180) / 40.0) ** 2)
    rss[0].data = level[:, numpy.newaxis] * shape
    # zero order: extra light in some fibers, marked in the mask
    mask = numpy.zeros(rss[0].data.shape, dtype="uint8")
    rss[0].data[100:300, 170:190] += 40.0
    mask[100:300, 170:190] = 1
    rss.append(flat_mask_hdu(mask))
    # a fiber invalid in the fiber flat, much brighter, is not used in the mean
    fibers = fits.ImageHDU(name="FIBERS")
    fibers.header["FIB010_V"] = False
    rss.append(fibers)
    rss[0].data[9] *= 3.0

    inputs = make_input()
    inputs.normalize_region = []
    # the window starts at the centre of the spectrum of the test, with the zero order
    recipe.default_normalize_region = (100, 200)
    master = recipe.run(inputs).master_twilightflat.open()
    good = numpy.ones(623, dtype=bool)
    good[[4, 9]] = False
    # the mean of each fiber in the window, without the zero order, is its level times the same factor
    assert numpy.allclose(master[0].data[good], (level / level[good].mean())[good, numpy.newaxis])
    assert "FLATMASK" not in master
    history = str(master[0].header["HISTORY"])
    start, end = (int(v) for v in history.split("columns ")[1].split()[0].split("-"))
    # the signal grows to the right of the initial window, the window moves to just after the zero order
    assert (start, end) == (190, 290)


def test_twilight_recipe_fiber_without_coverage(recipe, twilight_rss):
    """A fiber without coverage, not listed in the calibration, does not block the region nor enter the mean"""
    level, rss = twilight_rss
    rss[0].data[622] = 0.0
    wlmap = numpy.ones(rss[0].data.shape)
    wlmap[622] = 0.0
    rss.append(fits.ImageHDU(wlmap, name="WLMAP"))
    # zero order in the initial window: the window has to move
    from megaradrp.processing.fiberflat import flat_mask_hdu

    mask = numpy.zeros(rss[0].data.shape, dtype="uint8")
    mask[100:300, 150:160] = 1
    rss.append(flat_mask_hdu(mask))

    inputs = make_input()
    inputs.normalize_region = []
    recipe.default_normalize_region = (100, 200)
    master = recipe.run(inputs).master_twilightflat.open()
    good = numpy.ones(623, dtype=bool)
    good[[4, 622]] = False
    collapse = level + 2.0
    assert numpy.allclose(master[0].data[good], (collapse / collapse[good].mean())[good, numpy.newaxis])
    assert numpy.isclose(master[0].data[good, 0].mean(), 1.0)
    assert "columns 160-260" in str(master[0].header["HISTORY"])
