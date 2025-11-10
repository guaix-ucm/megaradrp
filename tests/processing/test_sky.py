import numpy

from megaradrp.processing.sky import subtract_sky_rss
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
