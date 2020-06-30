import pytest
import astropy.io.fits as fits

from megaradrp.tests.simpleobj import create_spec_header
from ..wavecalibration import header_add_barycentric_correction


def test_add_barycentric_missing1():
    hdr = create_spec_header()
    del hdr['RADEG']

    with pytest.raises(KeyError):
        header_add_barycentric_correction(hdr, key='b')


def test_add_barycentric_missing2():
    hdr = create_spec_header()
    del hdr['DATE-OBS']

    with pytest.raises(KeyError):
        header_add_barycentric_correction(hdr, key='b')


def test_add_barycentric_missing3():
    hdr = fits.Header()
    hdr['DATE-OBS'] = '2017-08-23T21:38:30.55'
    # GTC
    hdr['OBSGEO-X'] = 5327285.0921
    hdr['OBSGEO-Y'] = -1718777.1125
    hdr['OBSGEO-Z'] = 3051786.7327

    hdr['RADEG'] = 285.481037748898
    hdr['DECDEG'] = 42.4882140636786

    with pytest.raises(TypeError):
        header_add_barycentric_correction(hdr, key='b')
