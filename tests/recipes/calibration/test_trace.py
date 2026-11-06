"""The plots of the trace recipe are saved without pyplot"""

import matplotlib.pyplot as plt
import numpy
import pytest

from megaradrp.recipes.calibration.trace import init_traces, refine_boxes_from_array


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
