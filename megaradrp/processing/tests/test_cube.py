
import pytest

import astropy.wcs

from megaradrp.tests.simpleobj import create_spec_header, create_sky_header
from megaradrp.processing.wavecalibration import header_add_barycentric_correction

from ..cube import create_cube, merge_wcs, merge_wcs_alt


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)


def test_merge_wcs():
    hdr1 = create_spec_header()
    hdr1 = header_add_barycentric_correction(hdr1)
    hdr2 = create_sky_header()
    res = merge_wcs(hdr2, hdr1)
    cunit3 = res['CUNIT3']
    assert cunit3 == ''


def test_merge_wcs_2():
    import astropy.wcs
    hdr_sky = create_sky_header()
    hdr_spec = create_spec_header()
    hdr_spec = header_add_barycentric_correction(hdr_spec)
    allw = astropy.wcs.find_all_wcs(hdr_spec)
    out = hdr_spec.copy()
    for w in allw:
        ss = w.wcs.alt
        merge_wcs_alt(hdr_sky, hdr_spec, out, spec_suffix=ss)

    assert True


def test_merge2_wcs():
    hdr_sky = create_sky_header()
    hdr_spec = create_spec_header()
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