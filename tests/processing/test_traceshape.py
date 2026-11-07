import numpy
import pytest

from megaradrp.processing.traceshape import TraceShapeModel, complete_traces
from megaradrp.products.tracemap import GeometricTrace

XREF = 2000.0
XLIMITS = (6, 4090)
POLDEG = 5


def true_trace(y0, x):
    """A family of traces whose shape changes smoothly with the position y0"""
    u = (x - XREF) / 2000.0
    return y0 + (2.0 + 4.0 * y0 / 4000.0) * u + (8.0 - 10.0 * y0 / 4000.0) * u**2


def make_traces(nfibers=40, sep=6.0, first=1000.0, xmax=None):
    traces, points = [], {}
    for fibid in range(1, nfibers + 1):
        y0 = first + sep * (fibid - 1)
        x = numpy.arange(XLIMITS[0], XLIMITS[1] + 1, 2.0)
        if xmax is not None and fibid in xmax:
            x = x[x <= xmax[fibid]]
        y = true_trace(y0, x)
        pfit = numpy.polynomial.polynomial.polyfit(x, y, POLDEG)
        trace = GeometricTrace(fibid=fibid, boxid=0, start=int(x[0]), stop=int(x[-1]), fitparms=pfit.tolist())
        traces.append(trace)
        points[fibid] = (x, y)
    return traces, points


def test_shape_model_predicts_shape():
    xgrid = numpy.arange(6, 4091, 20.0)
    yref = numpy.linspace(100.0, 4000.0, 50)
    shapes = numpy.array([true_trace(y0, xgrid) - y0 for y0 in yref])
    model = TraceShapeModel(xgrid, yref, shapes, deg=4)
    for y0 in (150.0, 2222.0, 3950.0):
        expected = true_trace(y0, xgrid) - y0
        assert model.shape(y0, xgrid) == pytest.approx(expected, abs=1e-6)


def test_shape_model_needs_enough_traces():
    xgrid = numpy.arange(0, 100, 10.0)
    with pytest.raises(ValueError):
        TraceShapeModel(xgrid, [1.0, 2.0, 3.0], numpy.zeros((3, len(xgrid))), deg=4)


def test_complete_traces_truncated():
    """A trace that stops early is completed with the shape of the others"""
    traces, points = make_traces(xmax={25: 2600.0})
    before = traces[24]
    changed = complete_traces(traces, points, XLIMITS, XREF, POLDEG)
    assert changed == [25]
    after = traces[24]
    assert after is not before
    assert (after.start, after.stop) == XLIMITS
    x = numpy.arange(XLIMITS[0], XLIMITS[1] + 1, 10.0)
    assert after.polynomial(x) == pytest.approx(true_trace(1000.0 + 6.0 * 24, x), abs=0.05)


def test_complete_traces_jump():
    """The points of a trace that jumps to the next fiber are rejected"""
    traces, points = make_traces()
    x, y = points[30]
    y_jump = numpy.where(x > 3200, y + 6.0, y)
    points[30] = (x, y_jump)
    traces[29] = GeometricTrace(
        fibid=30,
        boxid=0,
        start=int(x[0]),
        stop=int(x[-1]),
        fitparms=numpy.polynomial.polynomial.polyfit(x, y_jump, POLDEG).tolist(),
    )
    changed = complete_traces(traces, points, XLIMITS, XREF, POLDEG)
    assert changed == [30]
    xx = numpy.arange(XLIMITS[0], XLIMITS[1] + 1, 10.0)
    assert traces[29].polynomial(xx) == pytest.approx(true_trace(1000.0 + 6.0 * 29, xx), abs=0.05)


def test_complete_traces_unchanged():
    """Complete traces are not changed"""
    traces, points = make_traces()
    before = list(traces)
    assert complete_traces(traces, points, XLIMITS, XREF, POLDEG) == []
    assert all(a is b for a, b in zip(traces, before))


def test_complete_traces_too_few():
    """Without enough complete traces, nothing is changed"""
    traces, points = make_traces(nfibers=10, xmax={5: 2600.0})
    assert complete_traces(traces, points, XLIMITS, XREF, POLDEG, min_complete=20) == []
    assert traces[4].stop == 2600
