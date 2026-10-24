import functools

import numpy
import pytest

import megaradrp.instrument.constants as cons
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.centroid import calc_centroid, calc_centroid_brightest, find_brightest_spaxel
from megaradrp.testing.create_image import create_rss

PLATESCALE = cons.GTC_FC_A_PLATESCALE.value


@functools.cache
def _fiber_header(unit):
    """Header of FIBERS of the LCB, with the positions in mm or in arcsec"""
    hdrf = create_default_fiber_header("LCB")
    if unit == "arcsec":
        for key in list(hdrf.keys()):
            if key.endswith("_X") or key.endswith("_Y"):
                hdrf[key] = hdrf[key] * PLATESCALE
        hdrf["FUNIT"] = "arcsec"
    return hdrf, FocalPlaneConf.from_header(hdrf)


def fiber_header(unit):
    return _fiber_header(unit)[0].copy()


def star_rss(hdrf, fibid, offset=(0.0, 0.0), sigma_mm=0.5):
    """RSS with a Gaussian star near the fiber fibid, offset in mm"""
    unit = hdrf.get("FUNIT", "mm")
    fp_conf = _fiber_header(unit)[1]
    scale = PLATESCALE if unit == "arcsec" else 1.0
    center = fp_conf.fibers[fibid]
    x0, y0 = center.x + offset[0] * scale, center.y + offset[1] * scale
    data = numpy.zeros((623, 100), dtype="float32")
    for fiber in fp_conf.connected_fibers():
        dist2 = ((fiber.x - x0) ** 2 + (fiber.y - y0) ** 2) / scale**2
        data[fiber.fibid - 1] = 100.0 * numpy.exp(-dist2 / (2 * sigma_mm**2))
    rss = create_rss(data, numpy.ones_like(data))
    rss["FIBERS"].header = hdrf
    return rss, center


def test_find_brightest_spaxel():
    rss, _ = star_rss(fiber_header("mm"), 311)
    fiber = find_brightest_spaxel(rss, (10, 90))
    assert fiber.fibid == 311


@pytest.mark.parametrize("unit, nrings", [("mm", 1), ("mm", 2), ("mm", 3), ("arcsec", 2)])
def test_calc_centroid(unit, nrings):
    """The centroid of a star centered in a fiber, in mm whatever the unit of the fibers"""
    rss, center = star_rss(fiber_header(unit), 311)
    centroid = calc_centroid(rss, (10, 90), (center.x, center.y), nrings)
    scale = PLATESCALE if unit == "arcsec" else 1.0
    assert numpy.allclose(centroid, (center.x / scale, center.y / scale))


def test_calc_centroid_shifted():
    """The centroid moves towards the star"""
    hdrf = fiber_header("mm")
    rss, center = star_rss(hdrf, 311, offset=(0.15, 0.0))
    centroid = calc_centroid(rss, (10, 90), (center.x, center.y), 2)
    assert center.x < centroid[0] < center.x + 0.15
    assert centroid[1] == pytest.approx(center.y, abs=0.01)


@pytest.mark.parametrize("unit", ["mm", "arcsec"])
def test_calc_centroid_brightest(unit):
    rss, center = star_rss(fiber_header(unit), 311)
    centroid = calc_centroid_brightest(rss, (10, 90), 2)
    scale = PLATESCALE if unit == "arcsec" else 1.0
    assert numpy.allclose(centroid, (center.x / scale, center.y / scale))
