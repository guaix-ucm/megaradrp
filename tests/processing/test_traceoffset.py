import math
import types

import numpy
import pytest
from numpy.polynomial import Polynomial

from megaradrp.processing.traceoffset import (
    TraceOffset,
    _base_sigma,
    _pixel_gaussian,
    combine_extraction_offsets,
    measure_extraction_offset,
)
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


def fake_result(c0, c1=0.0, reliable=True):
    return TraceOffset(
        coef=[c0, c1],
        error=[0.01, 0.0],
        reliable=reliable,
        alias_contrast=100.0,
        alias_bootstrap=1.0 if reliable else 0.6,
        nsignal=40,
        sigma_scale=1.3,
        reduced_chi2=1.0,
        bands=[],
        grid=numpy.zeros(1),
        curve=numpy.zeros(1),
        message="" if reliable else "only 2 fibers with signal",
    )


def test_combine_extraction_offsets(caplog):
    """The mean of the reliable measurements; the others are reported and left out"""
    tracemap = make_tracemap()
    results = [fake_result(1.0, 1e-4), fake_result(1.2, 3e-4), fake_result(5.0, reliable=False)]
    coef, history = combine_extraction_offsets(results, tracemap)
    assert coef == pytest.approx([1.1, 2e-4])
    assert len(history) == 4
    assert "not reliable (only 2 fibers with signal)" in history[2]
    assert history[3] == "Extraction offset measured: [1.1000, 0.0002]"
    assert "reliable in 2 of 3 images" in caplog.text
    assert "differ" not in caplog.text


def test_combine_extraction_offsets_spread(caplog):
    """Images with offsets more different than max_spread give a warning"""
    coef, history = combine_extraction_offsets([fake_result(1.0), fake_result(1.5)], make_tracemap(), max_spread=0.3)
    assert coef == pytest.approx([1.25, 0.0])
    assert "differ by 0.50 pixels" in caplog.text
    assert any("differ by 0.50 pixels" in line for line in history)


def test_combine_extraction_offsets_none_reliable(caplog):
    coef, history = combine_extraction_offsets([fake_result(3.0, reliable=False)], make_tracemap())
    assert coef == [0.0]
    assert history[-1] == "No reliable measurement of the extraction offset, using 0"
    assert "no reliable measurement" in caplog.text
