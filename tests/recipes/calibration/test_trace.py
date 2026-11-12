"""The plots of the trace recipe are saved without pyplot"""

import matplotlib.pyplot as plt
import numpy
import pytest

from megaradrp.recipes.calibration.trace import (
    box_contrast,
    flag_traces,
    init_traces,
    refine_boxes_from_array,
    select_reference_column,
)


def test_refine_boxes_from_array(tmp_path, monkeypatch):
    """The frontiers between the boxes are found 7 pixels after the expected ones"""
    monkeypatch.chdir(tmp_path)
    expected = [300, 600, 900, 1200, 1500]
    y = numpy.arange(1800.0)
    profile = 1000.0 + 50 * numpy.cos(y / 5)
    for frontier in expected:
        profile -= 800 * numpy.exp(-0.5 * ((y - frontier - 7) / 3) ** 2)
    arr = numpy.repeat(profile[:, numpy.newaxis], 20, axis=1)
    before = plt.get_fignums()
    refined, col = refine_boxes_from_array(arr, expected, 10, intermediate_results=True)
    assert refined == [frontier + 7 for frontier in expected]
    assert col == 10
    names = {path.name for path in tmp_path.iterdir()}
    assert names == {"offset_borders_cross.png", "frontiers_between_pseudoslits.png"}
    assert plt.get_fignums() == before


def _column_image(positions, nrows=200, ncols=31, sigma=1.5, level=10000.0):
    y = numpy.arange(nrows)[:, numpy.newaxis]
    profile = numpy.zeros((nrows, 1))
    for p in positions:
        profile = profile + level * numpy.exp(-0.5 * ((y - p) / sigma) ** 2)
    return numpy.repeat(profile + 100.0, ncols, axis=1)


BOXES = [{"name": "a", "nfibers": 10}, {"name": "b", "nfibers": 10}]
BORDERS = [10, 80, 146]
FIBERS_A = [16.0 + 6 * k for k in range(10)]
FIBERS_B = [86.0 + 6 * k for k in range(10)]


def _identified(traces):
    return {t.fibid: (None if t.start is None else round(t.start[1], 1)) for t in traces}


def test_init_traces_box_with_missing_fibers(tmp_path, monkeypatch):
    """A box without its last fibers does not shift the identification of the next box"""
    monkeypatch.chdir(tmp_path)
    image = _column_image(FIBERS_A[:7] + FIBERS_B)
    traces = init_traces(image, center=15, hs=3, boxes=BOXES, box_borders=BORDERS, threshold=0.3)
    found = _identified(traces)
    assert [found[f] for f in range(8, 11)] == [None, None, None]
    for fibid, pos in zip(range(1, 8), FIBERS_A):
        assert found[fibid] == pytest.approx(pos, abs=0.1)
    for fibid, pos in zip(range(11, 21), FIBERS_B):
        assert found[fibid] == pytest.approx(pos, abs=0.1)
    # without debug_plot, no figures are saved
    assert list(tmp_path.iterdir()) == []


def test_init_traces_box_without_fibers(tmp_path, monkeypatch):
    """A box without fibers keeps its fiber ids, and the next box is numbered after it"""
    monkeypatch.chdir(tmp_path)
    image = _column_image(FIBERS_B)
    traces = init_traces(image, center=15, hs=3, boxes=BOXES, box_borders=BORDERS, threshold=0.3)
    found = _identified(traces)
    assert sorted(found) == list(range(1, 21))
    assert all(found[f] is None for f in range(1, 11))
    for fibid, pos in zip(range(11, 21), FIBERS_B):
        assert found[fibid] == pytest.approx(pos, abs=0.1)


def test_init_traces_missing_fiber_inside_box(tmp_path, monkeypatch):
    """A dark fiber inside a box does not lose the rest of the box"""
    monkeypatch.chdir(tmp_path)
    image = _column_image(FIBERS_A[:3] + FIBERS_A[4:] + FIBERS_B)
    traces = init_traces(image, center=15, hs=3, boxes=BOXES, box_borders=BORDERS, threshold=0.3)
    found = _identified(traces)
    assert found[4] is None
    for fibid, pos in zip(range(1, 11), FIBERS_A):
        if fibid != 4:
            assert found[fibid] == pytest.approx(pos, abs=0.1)


def _defocused_image(ncols=4096, nrows=200):
    """Two boxes of fibers; the fibers of the second box widen towards the red"""
    y = numpy.arange(nrows)[:, numpy.newaxis]
    x = numpy.arange(ncols)[numpy.newaxis, :]
    image = numpy.full((nrows, ncols), 100.0)
    for pos in FIBERS_A:
        image += 10000.0 * numpy.exp(-0.5 * ((y - pos) / 1.5) ** 2)
    sigma = 1.5 + 2.5 * numpy.clip((x - 600.0) / 2000.0, 0, 1)
    for pos in FIBERS_B:
        image += 10000.0 * numpy.exp(-0.5 * ((y - pos) / sigma) ** 2)
    return image


def test_box_contrast():
    """The contrast of a box is lower for defocused fibers"""
    image = _defocused_image()
    blue = box_contrast(image[:, 400], BORDERS)
    red = box_contrast(image[:, 3500], BORDERS)
    assert blue[0] > 0.6
    assert blue[1] == pytest.approx(blue[0])
    assert red[0] == pytest.approx(blue[0])
    assert red[1] < 0.1


def test_select_reference_column():
    """The selected column is where the defocused box is still sharp"""
    image = _defocused_image()
    candidates = range(400, 3700, 100)
    best, scores = select_reference_column(image, BORDERS, candidates, hw=5)
    assert best <= 600
    assert len(scores) == len(candidates)
    assert max(scores) == scores[list(candidates).index(best)]


def _trace(fibid, y0, start=6, stop=4090, slope=0.0):
    from megaradrp.products.tracemap import GeometricTrace

    return GeometricTrace(fibid=fibid, boxid=0, start=start, stop=stop, fitparms=[y0, slope])


def test_flag_traces():
    """Weak, truncated and crossing traces are flagged; truncated ones are invalid after completion"""
    traces = [_trace(f, 100.0 + 6 * f) for f in range(1, 21)]
    # fiber 5 is weak
    peaks = {f: 1000.0 for f in range(1, 21)}
    peaks[5] = 100.0
    # fiber 10 stops at 2500
    traces[9] = _trace(10, 160.0, stop=2500)
    # fiber 15 crosses fiber 16
    traces[14] = _trace(15, 190.0, slope=0.01)
    flag_traces(traces, peaks, (6, 4090), completed=True)
    flags = {t.fibid: t.flags for t in traces}
    assert flags[5] == ["weak"]
    assert flags[10] == ["truncated", "not_completed"]
    assert flags[15] == ["crossing"]
    assert flags[16] == ["crossing"]
    assert sum(1 for f in flags.values() if f) == 4
    assert traces[4].valid
    assert not traces[9].valid
    assert not traces[14].valid


def test_flag_traces_without_completion():
    """Without completion, a truncated trace is only flagged"""
    traces = [_trace(f, 100.0 + 6 * f) for f in range(1, 21)]
    traces[9] = _trace(10, 160.0, stop=2500)
    flag_traces(traces, {f: 1000.0 for f in range(1, 21)}, (6, 4090), completed=False)
    assert traces[9].flags == ["truncated"]
    assert traces[9].valid


def test_search_traces_skips_non_finite_points(monkeypatch):
    """A point of the tracer with NaN position does not make the polynomial NaN"""
    import astropy.io.fits as fits

    import megaradrp.recipes.calibration.trace as trace_mod

    real_trace = trace_mod.trace_func

    def trace_with_nan(*args, **kwargs):
        # numina < 0.42 can end a trace with a point with NaN position
        mm = real_trace(*args, **kwargs)
        return numpy.vstack([mm, [mm[-1, 0] + 2, numpy.nan, numpy.nan]])

    monkeypatch.setattr(trace_mod, "trace_func", trace_with_nan)
    image = _column_image(FIBERS_A + FIBERS_B, ncols=400)
    reduced = fits.HDUList([fits.PrimaryHDU(image)])
    recipe = trace_mod.TraceMapRecipe()
    contents, error_fitting, missing = recipe.search_traces(
        reduced, BOXES, BORDERS, cstart=200, poldeg=3, xmax_fit=400, complete=False
    )
    assert error_fitting == [] and missing == []
    assert len(contents) == 20
    for trace, pos in zip(contents, FIBERS_A + FIBERS_B):
        assert trace.valid
        assert trace.polynomial(200) == pytest.approx(pos, abs=0.1)
