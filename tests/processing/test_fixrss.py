import math

from astropy.io import fits
import numpy as np
import pytest

import megaradrp.instrument.constants as cons
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.fixrss import fix_missing_fiber, recompute_wcs
from megaradrp.testing.create_header import create_sky_header2

MISSING = 623


@pytest.fixture(scope="module")
def lcb_fiber_header():
    return create_default_fiber_header("LCB")


@pytest.fixture(scope="module")
def nearby(lcb_fiber_header):
    return FocalPlaneConf.from_header(lcb_fiber_header).nearby_fibers(MISSING)


def make_rss(nearby, fiber_header, nwl=20):
    """RSS where the fibers around MISSING have values 1, 2, 3..."""
    data = np.zeros((623, nwl), dtype="float32")
    for value, fibid in enumerate(nearby, start=1):
        data[fibid - 1] = value
    return fits.HDUList([fits.PrimaryHDU(data), fits.ImageHDU(header=fiber_header.copy(), name="FIBERS")])


def test_fix_missing_fiber(lcb_fiber_header, nearby):
    """The missing fiber is the average of the fibers around it"""
    rss = make_rss(nearby, lcb_fiber_header)
    result = fix_missing_fiber(rss, MISSING)
    expected = np.mean(range(1, len(nearby) + 1))
    assert np.allclose(result[0].data[MISSING - 1], expected)


def test_fix_missing_fiber_limits(lcb_fiber_header, nearby):
    """The limits of the missing fiber are the narrowest of the fibers around it"""
    rss = make_rss(nearby, lcb_fiber_header)
    hdr = rss["FIBERS"].header
    for idx, fibid in enumerate(nearby):
        hdr[f"FIB{fibid:03d}W1"] = 3 + idx
        hdr[f"FIB{fibid:03d}W2"] = 15 - idx
    w1 = 3 + len(nearby) - 1
    w2 = 15 - len(nearby) + 1

    result = fix_missing_fiber(rss, MISSING)
    assert hdr[f"FIB{MISSING:03d}W1"] == w1
    assert hdr[f"FIB{MISSING:03d}W2"] == w2
    row = result[0].data[MISSING - 1]
    # the pixels outside the limits are zero
    assert np.all(row[: w1 - 1] == 0)
    assert np.all(row[w2 - 1 :] == 0)
    assert np.all(row[w1 - 1 : w2 - 1] > 0)


@pytest.mark.parametrize("ipa", [0.0, 30.0, 200.0])
def test_recompute_wcs(ipa):
    hdr = create_sky_header2()
    hdr["CDELT1"] = abs(hdr["CDELT1"])
    result = recompute_wcs(hdr, ipa)
    pa = math.radians(ipa - cons.MEGARA_IAA.value)
    assert result["PC1_1"] == pytest.approx(math.cos(pa))
    assert result["PC2_2"] == pytest.approx(math.cos(pa))
    assert result["PC1_2"] == pytest.approx(-math.sin(pa))
    assert result["PC2_1"] == pytest.approx(math.sin(pa))
    # the axis of RA increases to the left
    assert result["CDELT1"] < 0
