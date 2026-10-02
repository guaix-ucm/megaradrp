import warnings

import astropy.wcs
import numpy
import pytest

from megaradrp.testing.create_header import create_spec_header, create_spec_header2
from megaradrp.processing.wavecalibration import header_add_barycentric_correction

OBSGEO_XYZ = ["OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"]
OBSGEO_BLH = ["OBSGEO-L", "OBSGEO-B", "OBSGEO-H"]


def check_barycentric(hdr):
    assert numpy.isfinite(hdr["VELOSYSB"])
    assert numpy.isfinite(hdr["CRVAL1B"])
    assert numpy.isfinite(hdr["CDELT1B"])
    # The correction is at most ~30 km/s
    assert abs(hdr["VELOSYSB"]) < 35000
    assert hdr["SPECSYSB"] == "BARYCENT"


def header_without(keys):
    hdr = create_spec_header2()
    for key in keys:
        del hdr[key]
    return hdr


def reference_velocity():
    return header_add_barycentric_correction(create_spec_header2())["VELOSYSB"]


def test_barycentric_xyz_and_blh():
    hdr = header_add_barycentric_correction(create_spec_header2())
    check_barycentric(hdr)


def test_barycentric_xyz_only():
    hdr = header_add_barycentric_correction(header_without(OBSGEO_BLH))
    check_barycentric(hdr)
    assert hdr["VELOSYSB"] == pytest.approx(reference_velocity(), abs=1.0)


def test_barycentric_blh_only():
    """OBSGEO-B/L/H only, as written by previous versions: no default values"""
    hdr = header_add_barycentric_correction(header_without(OBSGEO_XYZ))
    check_barycentric(hdr)
    assert "OBSGEO-X" not in hdr
    # less than 1 m/s
    assert hdr["VELOSYSB"] == pytest.approx(reference_velocity(), abs=1.0)
    # the header used in other tests
    check_barycentric(header_add_barycentric_correction(create_spec_header()))


def test_barycentric_no_obsgeo():
    """Without OBSGEO keywords, the geocentric coordinates of GTC are added"""
    with pytest.warns(RuntimeWarning, match="OBSGEO- keywords not defined"):
        hdr = header_add_barycentric_correction(header_without(OBSGEO_XYZ + OBSGEO_BLH))
    check_barycentric(hdr)
    assert hdr["OBSGEO-X"] == 5327285.0921
    for key in OBSGEO_BLH:
        assert key not in hdr
    assert hdr["VELOSYSB"] == pytest.approx(reference_velocity(), abs=1.0)


def test_default_obsgeo_consistent():
    """wcslib does not find the default coordinates inconsistent"""
    with pytest.warns(RuntimeWarning, match="OBSGEO- keywords not defined"):
        hdr = header_add_barycentric_correction(header_without(OBSGEO_XYZ + OBSGEO_BLH))
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        astropy.wcs.WCS(hdr, fix=True)
    messages = [str(w.message) for w in record]
    assert not [m for m in messages if "inconsistent" in m], messages
