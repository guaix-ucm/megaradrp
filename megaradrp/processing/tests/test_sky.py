import numpy

from megaradrp.datamodel import create_default_fiber_header
from ..sky import subtract_sky_rss


def create_rss(value, wlmap):
    import astropy.io.fits as fits
    data1 = value + numpy.zeros((623, 4300), dtype='float32')
    hdu = fits.PrimaryHDU(data1)
    hdrf = create_default_fiber_header('LCB')
    fibers = fits.ImageHDU(header=hdrf, name='FIBERS')
    rss_map = fits.ImageHDU(wlmap, name='WLMAP')
    return fits.HDUList([hdu, fibers, rss_map])


def test_subtract_sky_rss():

    wlmap = numpy.zeros((623, 4300), dtype='float32')
    wlmap[:,350:4105] = 1.0
    wlmap[622,:] = 0
    img1 = create_rss(1000, wlmap)
    img2 = create_rss(400, wlmap)

    final_img, img, sky_img = subtract_sky_rss(img1, img2)
    assert img is img1
    # In final image, regions outside WLMAP must be at zero
    assert final_img[0].data[:, 100:200].min() == 0
    assert final_img[0].data[:, 100:200].max() == 0

    assert final_img[0].data[622, :].max() == 0
    assert final_img[0].data[622, :].min() == 0
