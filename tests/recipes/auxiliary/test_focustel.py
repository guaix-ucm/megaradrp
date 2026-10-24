import types

from astropy.io import fits
import numpy
import pytest
from numina.core import DataFrame

import megaradrp.recipes.auxiliary.focustel as focustel
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.testing.create_image import create_rss

BEST_FOCUS = 20.0


@pytest.fixture(scope="module")
def lcb():
    hdrf = create_default_fiber_header("LCB")
    return hdrf, FocalPlaneConf.from_header(hdrf)


def star_rss(lcb, focus):
    """RSS of a star in the center of the LCB, wider away from the best focus"""
    hdrf, fp_conf = lcb
    sigma = 0.4 + 0.1 * (focus - BEST_FOCUS) ** 2
    data = numpy.zeros((623, 3200), dtype="float32")
    for fiber in fp_conf.connected_fibers():
        data[fiber.fibid - 1] = 100.0 * numpy.exp(-(fiber.x**2 + fiber.y**2) / (2 * sigma**2))
    rss = create_rss(data, numpy.ones_like(data, dtype="int16"))
    rss["FIBERS"].header = hdrf.copy()
    rss[0].header["M2UZ"] = focus
    return rss


@pytest.fixture
def recipe(monkeypatch, lcb):
    """The recipe, without the 2D processing, the extraction and the intermediate files"""
    recipe = focustel.FocusTelescopeRecipe()
    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: None)
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: None)
    # the combination of a group of frames is its first frame
    monkeypatch.setattr(
        focustel, "basic_processing_with_combination_frames", lambda frames, flow, **kwds: frames[0].open()
    )
    # the extracted RSS has a star whose width depends on the focus
    monkeypatch.setattr(
        recipe, "run_reduction_1d", lambda img, *args, **kwds: (img, star_rss(lcb, img[0].header["M2UZ"]))
    )
    return recipe


def make_input(focus_values, **kwds):
    frames = []
    for focus in focus_values:
        hdu = fits.PrimaryHDU()
        hdu.header["M2UZ"] = focus
        frames.append(DataFrame(frame=fits.HDUList([hdu])))
    values = dict(
        obresult=types.SimpleNamespace(frames=frames, configuration=None),
        position=(0, 0),
        master_apertures=None,
        master_wlcalib=None,
        master_fiberflat=None,
        extraction_offset=[0.0],
        sky_rss=None,
        ignored_sky_bundles=[],
    )
    values.update(kwds)
    return types.SimpleNamespace(**values)


@pytest.mark.parametrize("ignored_sky_bundles", [[], pytest.param([93, 94], marks=pytest.mark.slow)])
def test_focus_telescope(recipe, ignored_sky_bundles):
    """The best focus, with the sky subtracted from the SKY bundles (#352)"""
    # three groups of images, the minimum for the fit
    focus_values = [19.0, 20.0, 21.0, 19.0]
    result = recipe.run(make_input(focus_values, ignored_sky_bundles=ignored_sky_bundles))
    assert result.focus_table == pytest.approx(BEST_FOCUS, abs=1e-6)


def test_focus_telescope_one_focus(recipe):
    from numina.exceptions import RecipeError

    with pytest.raises(RecipeError):
        recipe.run(make_input([20.0, 20.0]))
