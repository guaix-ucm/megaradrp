import math

import pytest

from ..cube import create_cube, merge_wcs


def generate_wcs1():
    import astropy.wcs

    wcsl = astropy.wcs.WCS(naxis=2)

    crpix = 3.0
    crval = 5300.0
    cdelt = 0.1

    wcsl.wcs.crpix = [crpix, 0.0]
    wcsl.wcs.crval = [crval, 0.0]
    wcsl.wcs.cdelt = [cdelt, 1.0]
    wcsl.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    ang = math.pi / 3.0
    wcsl.wcs.pc = [[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]]
    return wcsl


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)


def test_merge_wcs():
    hdr1 = generate_wcs1().to_header()
    hdr2 = generate_wcs1().to_header()
    res = merge_wcs(hdr1, hdr2)
    cunit3 = res['CUNIT3']
    assert cunit3 == ''