import math
import types

import numpy
import pytest
from numpy.polynomial import Polynomial

from megaradrp.processing.traceoffset import _base_sigma, _pixel_gaussian, measure_extraction_offset
from megaradrp.products.tracemap import GeometricTrace, TraceMap
from megaradrp.testing.create_traceoffset import SIGMA
from megaradrp.testing.create_traceoffset import create_offset_image as make_image
from megaradrp.testing.create_traceoffset import create_offset_tracemap as make_tracemap

OPTS = dict(nbands=4, band_width=32, search=(-8.0, 8.0), step=0.25)


def test_pixel_gaussian():
    rows = numpy.arange(-20, 21)
    values = _pixel_gaussian(rows, 0.3, 1.5)
    assert math.isclose(values.sum(), 1.0, rel_tol=1e-9)
    assert numpy.argmax(values) == 20


@pytest.mark.parametrize("c0", [2.3, -1.1, -5.5])
def test_constant_offset(c0):
    """Also an offset close to the distance between fibers: the gaps between blocks break the alias"""
    tracemap = make_tracemap()
    result = measure_extraction_offset(make_image(tracemap, c0), tracemap, fit_slope=False, **OPTS)
    assert result.reliable, result.message
    assert result.coef[0] == pytest.approx(c0, abs=0.03)
    assert result.coef[1] == 0.0
    assert result.nsignal == 40
    assert result.sigma_scale == pytest.approx(SIGMA, abs=0.25)
    assert all(band["c0"] == pytest.approx(c0, abs=0.05) for band in result.bands)


def test_linear_offset():
    tracemap = make_tracemap()
    c0, c1 = 1.0, 0.003
    result = measure_extraction_offset(make_image(tracemap, c0, c1), tracemap, **OPTS)
    assert result.reliable, result.message
    y = numpy.array([30.0, 290.0])
    assert numpy.allclose(result.coef[0] + result.coef[1] * y, c0 + c1 * y, atol=0.05)


def test_global_offset_applied():
    """The offset is measured on top of the global_offset of the calibration"""
    tracemap = make_tracemap()
    image = make_image(tracemap, 2.0)
    tracemap.global_offset = Polynomial([0.5])
    result = measure_extraction_offset(image, tracemap, fit_slope=False, **OPTS)
    assert result.coef[0] == pytest.approx(1.5, abs=0.03)


def test_few_fibers_no_slope():
    """Fibers with signal only at one end of the slit: no slope"""
    tracemap = make_tracemap()
    flux = numpy.zeros(40)
    flux[:8] = 3000.0
    result = measure_extraction_offset(make_image(tracemap, 1.7, flux=flux), tracemap, **OPTS)
    assert result.coef[1] == 0.0
    assert result.coef[0] == pytest.approx(1.7, abs=0.05)
    assert result.nsignal >= 8


def test_no_signal():
    tracemap = make_tracemap()
    result = measure_extraction_offset(make_image(tracemap, 1.0, flux=numpy.zeros(40)), tracemap, **OPTS)
    assert not result.reliable
    assert result.nsignal < 3
    assert "fibers with signal" in result.message


def test_mask():
    """Masked pixels, here a saturated column, are left out"""
    tracemap = make_tracemap()
    image = make_image(tracemap, 0.8)
    image[:, 100:104] = 1e6
    mask = numpy.zeros(image.shape, dtype=bool)
    mask[:, 100:104] = True
    result = measure_extraction_offset(image, tracemap, mask=mask, fit_slope=False, **OPTS)
    assert result.coef[0] == pytest.approx(0.8, abs=0.03)


def test_no_valid_traces():
    tracemap = TraceMap(instrument="MEGARA")
    tracemap.contents.append(GeometricTrace(1, 1, 0, 10, fitparms=[]))
    with pytest.raises(ValueError, match="no valid traces"):
        measure_extraction_offset(numpy.zeros((10, 10)), tracemap)


def test_base_sigma():
    """The width of the profiles of a ModelMap; 1 for a TraceMap"""
    model = types.SimpleNamespace(model={"params": {"stddev": lambda x: 1.0 + x / 100.0}})
    trace = types.SimpleNamespace()
    assert numpy.allclose(_base_sigma([model, trace], [0.0, 100.0]), [1.5, 1.0])
