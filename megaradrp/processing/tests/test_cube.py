
import pytest

from megaradrp.tests.simpleobj import create_spec_header, create_sky_header
from megaradrp.processing.wavecalibration import header_add_barycentric_correction

from ..cube import create_cube, merge_wcs


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)


def test_merge_wcs():
    hdr1 = create_sky_header()
    hdr2 = create_spec_header()
    hdr2 = header_add_barycentric_correction(hdr2)
    res = merge_wcs(hdr1, hdr2)
    cunit3 = res['CUNIT3']
    assert cunit3 == ''
