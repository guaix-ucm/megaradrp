import numpy
import pytest

from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.sky import subtract_sky, subtract_sky_rss
from megaradrp.testing.create_image import create_rss, create_scene_1212


def test_subtract_sky_rss():

    wlmap = numpy.zeros((623, 4300), dtype="float32")
    wlmap[:, 350:4105] = 1.0
    wlmap[622, :] = 0
    scene1 = create_scene_1212(1000)
    scene2 = create_scene_1212(400)
    img1 = create_rss(scene1, wlmap)
    img2 = create_rss(scene2, wlmap)

    final_img, img, sky_img = subtract_sky_rss(img1, img2)
    assert img is img1
    # In final image, regions outside WLMAP must be at zero
    assert final_img[0].data[:, 100:200].min() == 0
    assert final_img[0].data[:, 100:200].max() == 0

    assert final_img[0].data[622, :].max() == 0
    assert final_img[0].data[622, :].min() == 0
    # inside WLMAP, the sky image is subtracted
    assert numpy.allclose(final_img[0].data[:622, 350:4105], 600)


SKY_LEVEL = 5.0


@pytest.fixture(scope="module")
def focal_plane():
    return FocalPlaneConf.from_header(create_default_fiber_header("LCB"))


def make_rss(nwl=60):
    """RSS of the LCB with a constant sky, all the fibers covered in [10, 50)"""
    data = numpy.zeros((623, nwl), dtype="float32")
    wlmap = numpy.zeros((623, nwl), dtype="int16")
    wlmap[:, 10:50] = 1
    data[wlmap > 0] = SKY_LEVEL
    return create_rss(data, wlmap)


def test_subtract_sky(focal_plane):
    """The average of the sky fibers is subtracted from the valid fibers"""
    rss = make_rss()
    # an object in an IFU fiber
    rss[0].data[99, 10:50] += 100.0
    final_img, img, sky_img = subtract_sky(rss)

    assert img is rss
    assert numpy.allclose(final_img[0].data[99, 10:50], 100.0)
    ifu = [fiber.fibid - 1 for fiber in focal_plane.connected_fibers() if fiber.fibid != 100]
    assert numpy.allclose(final_img[0].data[ifu, 10:50], 0.0)
    sky_rows = [fibid - 1 for fibid in focal_plane.sky_fibers()]
    assert numpy.allclose(final_img[0].data[sky_rows, 10:50], 0.0)
    # outside the coverage of all the fibers, nothing is subtracted
    assert numpy.all(final_img[0].data[:, :10] == 0)
    # the sky image has the data of the sky fibers
    assert numpy.allclose(sky_img[0].data[sky_rows, 10:50], SKY_LEVEL)
    assert numpy.all(sky_img[0].data[99] == 0)
    # the input is not modified
    assert numpy.allclose(rss[0].data[99, 10:50], 100.0 + SKY_LEVEL)


def test_subtract_sky_ignored_bundles(focal_plane):
    """The fibers of the ignored sky bundles are not used"""
    rss = make_rss()
    bundle = 93
    for fibid in focal_plane.bundles[bundle].fibers:
        rss[0].data[fibid - 1, 10:50] += 70.0

    final_img, _, _ = subtract_sky(rss, ignored_sky_bundles=[bundle])
    assert numpy.allclose(final_img[0].data[0, 10:50], 0.0)

    # without ignoring it, the sky is biased
    final_img, _, _ = subtract_sky(rss)
    assert numpy.all(final_img[0].data[0, 10:50] < 0)


def test_subtract_sky_invalid_fibers():
    """The invalid fibers are not used as sky, and are not subtracted"""
    rss = make_rss()
    hdr = rss["FIBERS"].header
    fp_conf = FocalPlaneConf.from_header(hdr)
    sky_fibid = fp_conf.sky_fibers()[0]
    hdr[f"FIB{sky_fibid:03d}_V"] = False
    rss[0].data[sky_fibid - 1, 10:50] = 1000.0
    hdr["FIB100_V"] = False

    final_img, _, _ = subtract_sky(rss)
    assert numpy.allclose(final_img[0].data[0, 10:50], 0.0)
    assert numpy.allclose(final_img[0].data[99, 10:50], SKY_LEVEL)
