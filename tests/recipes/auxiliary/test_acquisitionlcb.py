import types

from astropy.io import fits
import numpy
import pytest

import megaradrp.recipes.scientific.base as base
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.recipes.auxiliary.acquisitionlcb import AcquireLCBRecipe
from megaradrp.testing.create_image import create_rss

NWL = 4000
SKY = 50.0


@pytest.fixture(scope="module")
def lcb():
    hdrf = create_default_fiber_header("LCB")
    return hdrf, FocalPlaneConf.from_header(hdrf)


def star_rss(lcb, x0, y0, sigma=0.4):
    """RSS with a star in (x0, y0), in mm, and a constant sky"""
    hdrf, fp_conf = lcb
    data = numpy.full((623, NWL), SKY, dtype="float32")
    for fiber in fp_conf.connected_fibers():
        dist2 = (fiber.x - x0) ** 2 + (fiber.y - y0) ** 2
        data[fiber.fibid - 1] += 1000 * numpy.exp(-dist2 / (2 * sigma**2))
    rss = create_rss(data, numpy.ones_like(data, dtype="int16"))
    rss["FIBERS"].header = hdrf.copy()
    rss[0].header["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    return rss


@pytest.fixture
def run_recipe(monkeypatch, lcb):
    """Run the recipe with the 2D reduction and the extraction replaced by a RSS"""

    def run(rss, **kwds):
        reduced2d = fits.HDUList([fits.PrimaryHDU(numpy.zeros((10, 10)))])
        monkeypatch.setattr(base.ImageRecipe, "base_run", lambda self, rinput: (reduced2d, rss))
        values = dict(ignored_sky_bundles=[], sky_rss=None, extraction_region=[1000, 3000], nrings=3, points=[(0, 0)])
        values.update(kwds)
        return AcquireLCBRecipe().run(types.SimpleNamespace(**values))

    return run


def test_acquisition_lcb(run_recipe, lcb):
    """The offset moves the star to the center of the LCB"""
    _, fp_conf = lcb
    # a star between spaxels, near the center
    center = min(fp_conf.connected_fibers(), key=lambda f: f.x**2 + f.y**2)
    x0, y0 = center.x + 0.2, center.y - 0.15
    result = run_recipe(star_rss(lcb, x0, y0))
    assert numpy.allclose(result.offset, [-x0, -y0], atol=0.05)
    # the sky is subtracted
    final = result.final_rss.open()[0].data
    reduced = result.reduced_rss.open()[0].data
    assert numpy.allclose(final[center.fibid - 1], reduced[center.fibid - 1] - SKY)
    assert result.reduced_image.open()[0].data.shape == (10, 10)


def test_acquisition_lcb_brightest(run_recipe, lcb):
    """The centroid is computed around the brightest spaxel"""
    _, fp_conf = lcb
    # a star centered in a spaxel, with the three rings inside the LCB
    fiber = min(fp_conf.connected_fibers(), key=lambda f: (f.x - 1.5) ** 2 + (f.y - 1.0) ** 2)
    rss = star_rss(lcb, fiber.x, fiber.y)
    # a fainter star far away
    other = min(fp_conf.connected_fibers(), key=lambda f: (f.x + 4) ** 2 + (f.y + 3) ** 2)
    rss[0].data[other.fibid - 1] += 200
    result = run_recipe(rss)
    assert numpy.allclose(result.offset, [-fiber.x, -fiber.y], atol=0.01)
