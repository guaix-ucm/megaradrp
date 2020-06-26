
import math

import astropy.wcs
import astropy.io.fits as fits


def create_simple_hdul():
    """Create a simple image for testing"""
    prim = fits.PrimaryHDU()
    prim.header['instrume'] = 'MEGARA'
    prim.header['VPH'] = 'LR-B'
    prim.header['DATE-OBS'] = '2019-02-21T01:02:02.2'
    prim.header['insmode'] = 'LCB'
    prim.header['UUID'] = '410f6ea0-c3df-481c-9820-5cf9e0ed3d91'
    fibers = fits.ImageHDU(name='FIBERS')
    fibers.header['CONFID'] = 'a908e9d2-7599-41f3-9e9d-91042df01da5'
    simple_img = fits.HDUList([prim, fibers])
    return simple_img


def generate_sky_wcs():
    wcsl = astropy.wcs.WCS(naxis=2)

    wcsl.wcs.crpix = [512, 512]
    wcsl.wcs.crval = [9.0000, 32.0000]
    wcsl.wcs.cdelt = [0.01, 0.01]
    wcsl.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    ang = math.pi / 3.0
    wcsl.wcs.pc = [[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]]
    return wcsl


def create_sky_header():
    wcsl = generate_sky_wcs()
    return wcsl.to_header()


def create_spec_header():
    hdr = fits.Header()
    hdr['DATE-OBS'] = '2017-08-23T21:38:30.55'
    # GTC
    hdr['OBSGEO-X'] = 5327285.0921
    hdr['OBSGEO-Y'] = -1718777.1125
    hdr['OBSGEO-Z'] = 3051786.7327

    hdr['RADEG'] = 285.481037748898
    hdr['DECDEG'] = 42.4882140636786

    hdr['CTYPE1'] = 'AWAV'
    hdr['CRPIX1'] = 1
    hdr['CRVAL1'] = 362.0
    hdr['CDELT1'] = 1.86
    hdr['CUNIT1'] = 'nm'

    hdr['CRPIX2'] = 0
    hdr['CRVAL2'] = 0
    hdr['CDELT2'] = 1
    hdr['CTYPE2'] = ''
    return hdr