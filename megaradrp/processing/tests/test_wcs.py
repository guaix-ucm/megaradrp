

import math

import numpy
import pytest

from ..wcs import compute_pa_from_ipa
from ..wcs import update_wcs_from_ipa


@pytest.mark.parametrize("ipa, pa", [(0, 163.854), (196.146,0), (90, 253.854)])
def test_pa_from_ipa1(ipa, pa):

    res = compute_pa_from_ipa(ipa)
    assert numpy.isclose([res % 360], [pa])


@pytest.mark.parametrize(
    "ipa, iaa, pa", [
        (0, -163.854, 163.854), (196.146, -163.854, 0), (90, -163.854, 253.854)
    ]
)
def test_pa_from_ipa2(ipa, iaa, pa):

    res = compute_pa_from_ipa(ipa, iaa)
    assert numpy.isclose([res % 360], [pa])


@pytest.mark.parametrize(
    "cdelt1", [
        0.1, -0.1
    ]
)
@pytest.mark.parametrize(
    "pa", [0.0, 45.0, 120.0,
    ]
)
def test_update_wcs_from_ipa(cdelt1, pa):

    hdr = {}
    hdr['CDELT1'] = cdelt1
    hdr = update_wcs_from_ipa(hdr, pa)
    pa_rad = numpy.deg2rad(pa)
    cos_a = math.cos(pa_rad)
    sin_a = math.sin(pa_rad)
    res = [[cos_a, -sin_a], [ sin_a, cos_a]]
    val = [[hdr['PC1_1'], hdr['PC1_2']], [hdr['PC2_1'], hdr['PC2_2']]]

    assert hdr['CDELT1'] == -abs(cdelt1)
    assert numpy.allclose(val, res)

