
import pytest
import astropy.io.fits as fits
import numpy

from ..multirss import generate_multi_rss

def generate_imgs(nimages, nfibers, nsamples):
    imgs = []
    for _ in range(nimages):
        data = numpy.empty((nfibers, nsamples))
        hdu1 = fits.PrimaryHDU(data)
        hdu2 = fits.ImageHDU(name='FIBERS')
        img = fits.HDUList([hdu1, hdu2])
        imgs.append(img)
    return imgs


def test_multirss():

    nimages = 10
    nfibers = 12
    nsamples = 300
    imgs = generate_imgs(nimages, nfibers, nsamples)

    result = generate_multi_rss(imgs)

    assert len(result) == nimages + 1

    primary = result[0]
    data = primary.data
    assert data.shape == (nimages * nfibers, nsamples)
    assert primary.header['MEG-NRSS'] == nimages
    # Check
    for idx, ext in enumerate(result[1:], 1):
        assert ext.header['EXTNAME'] == f"FIBERS{idx}"


def test_multirss_error():

    nimages = 0
    nfibers = 12
    nsamples = 300
    imgs = generate_imgs(nimages, nfibers, nsamples)

    with pytest.raises(ValueError):
        generate_multi_rss(imgs)
