
import pytest
import numpy as np
import astropy.wcs

from megaradrp.tests.simpleobj import create_spec_header2, create_sky_header2
from megaradrp.processing.wavecalibration import header_add_barycentric_correction

from ..cube import create_cube, merge_wcs


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)


def test_sub_wcs():
    hdr_sky = create_sky_header2()
    hdr_spec = create_spec_header2()
    hdr_spec = header_add_barycentric_correction(hdr_spec)
    wcs_sky = astropy.wcs.WCS(header=hdr_sky)
    wcs_spec = astropy.wcs.WCS(header=hdr_spec, key='B')
    wcs3 = wcs_sky.sub([1,2,0])
    wcs3.wcs.ctype[2] = wcs_spec.wcs.ctype[0]
    wcs3.wcs.crval[2] = wcs_spec.wcs.crval[0]
    wcs3.wcs.crpix[2] = wcs_spec.wcs.crpix[0]
    wcs3.wcs.cdelt[2] = wcs_spec.wcs.cdelt[0]
    wcs3.wcs.cunit[2] = wcs_spec.wcs.cunit[0]
    wcs3.wcs.specsys = wcs_spec.wcs.specsys
    wcs3.wcs.ssysobs = wcs_spec.wcs.ssysobs
    wcs3.wcs.velosys = wcs_spec.wcs.velosys
    hdr3 = wcs3.to_header(key='B')
    assert True


def test_merge_wcs():
    hdr_spec = create_spec_header2()
    hdr_spec = header_add_barycentric_correction(hdr_spec)
    hdr_sky = create_sky_header2()
    out = hdr_spec.copy()
    merge_wcs(hdr_sky, hdr_spec, out=out)

    wcs0 = astropy.wcs.WCS(header=out)
    wcsB = astropy.wcs.WCS(header=out, key='B')

    assert list(wcs0.wcs.ctype) == ['RA---TAN', 'DEC--TAN', 'AWAV']
    assert np.allclose([out['CRVAL1'], out['CRVAL2'], out['CRVAL3']],
                       [hdr_sky['CRVAL1'], hdr_sky['CRVAL2'], hdr_spec['CRVAL1']])
    assert np.allclose([out['CDELT1'], out['CDELT2'], out['CDELT3']],
                       [hdr_sky['CDELT1'], hdr_sky['CDELT2'], hdr_spec['CDELT1']])
    assert list(wcsB.wcs.ctype) == ['RA---TAN', 'DEC--TAN', 'AWAV']
    assert np.allclose([out['CRVAL1B'], out['CRVAL2B'], out['CRVAL3B']],
                       [hdr_sky['CRVAL1'], hdr_sky['CRVAL2'], hdr_spec['CRVAL1B']])
    assert np.allclose([out['CDELT1B'], out['CDELT2B'], out['CDELT3B']],
                       [hdr_sky['CDELT1'], hdr_sky['CDELT2'], hdr_spec['CDELT1B']])
