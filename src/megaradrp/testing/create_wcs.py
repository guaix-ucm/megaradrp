import math

import astropy.wcs


def generate_sky_wcs():
    wcsl = astropy.wcs.WCS(naxis=2)

    wcsl.wcs.crpix = [512, 512]
    wcsl.wcs.crval = [9.0000, 32.0000]
    wcsl.wcs.cdelt = [0.01, 0.01]
    wcsl.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    ang = math.pi / 3.0
    wcsl.wcs.pc = [[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]]
    return wcsl


def generate_wcs():
    wcsl = astropy.wcs.WCS(naxis=2)

    crpix = 3.0
    crval = 5300.0
    cdelt = 0.1

    wcsl.wcs.crpix = [crpix, 0.0]
    wcsl.wcs.crval = [crval, 0.0]
    wcsl.wcs.cdelt = [cdelt, 1.0]
    return wcsl
