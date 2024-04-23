

import pytest

from megaradrp.processing.fibermatch import generate_box_model, FiberModelElement
from megaradrp.processing.fibermatch import count_peaks, PeakMatch, PeakFound, PeakMode


PEAKS = [
    3.806000000000000000e+03,
    3.812000000000000000e+03,
    3.818000000000000000e+03,
    3.824000000000000000e+03,
    3.830000000000000000e+03,
    3.836000000000000000e+03,
    3.842000000000000000e+03,
    3.848000000000000000e+03,
    3.854000000000000000e+03,
    3.860000000000000000e+03,
    3.867000000000000000e+03,
    3.872000000000000000e+03,
    3.878000000000000000e+03,
    3.884000000000000000e+03,
    3.890000000000000000e+03,
    3.897000000000000000e+03,
    3.903000000000000000e+03,
    3.909000000000000000e+03,
    3.915000000000000000e+03,
    3.921000000000000000e+03
]


def test_generate_model():

    expected = [
        (1, PeakMode.FIBER_PEAK),
        (2, PeakMode.FIBER_PEAK),
        (3, PeakMode.FIBER_PEAK),
        (4, PeakMode.FIBER_PEAK),
        (5, PeakMode.FIBER_PEAK)
    ]
    model = generate_box_model(5, start=1)

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == FiberModelElement(fibid=e[0], mode=e[1])

    expected = [
        (1, PeakMode.FIBER_PEAK),
        (2, PeakMode.FIBER_DEAD),
        (3, PeakMode.FIBER_PEAK),
        (4, PeakMode.FIBER_PEAK),
        (5, PeakMode.FIBER_PEAK)
    ]
    model = generate_box_model(5, missing_relids=[2])

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == FiberModelElement(fibid=e[0], mode=e[1])

    expected = [
        (10, PeakMode.FIBER_PEAK),
        (12, PeakMode.FIBER_DEAD),
        (13, PeakMode.FIBER_PEAK),
        (14, PeakMode.FIBER_PEAK),
        (15, PeakMode.FIBER_PEAK)
    ]
    model = generate_box_model(5, start=10, skip_fibids=[
                               11], missing_relids=[2])

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == FiberModelElement(fibid=e[0], mode=e[1])


def test_count_peaks1():
    with pytest.raises(ValueError):
        count_peaks([])


def test_count_peaks():

    expected = []
    idx = 0
    for p in PEAKS:
        t = PeakMatch(idx + 1, p, PeakFound.FOUND_PEAK, idx)
        expected.append(t)
        idx += 1

    result = count_peaks(PEAKS, tol=1.2, distance=6.0)

    assert result == expected
