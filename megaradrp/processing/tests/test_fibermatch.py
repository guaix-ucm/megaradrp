

import pytest

from megaradrp.processing.fibermatch import generate_box_model
from megaradrp.processing.fibermatch import count_peaks


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
                   (1, 0),
                   (2, 0),
                   (3, 0),
                   (4, 0),
                   (5, 0)
               ]
    model = generate_box_model(5, start=1)

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == e

    expected = [
                   (1, 0),
                   (2, 1),
                   (3, 0),
                   (4, 0),
                   (5, 0)
               ]
    model = generate_box_model(5, missing_relids=[2])

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == e

    expected = [
                   (10, 0),
                   (12, 1),
                   (13, 0),
                   (14, 0),
                   (15, 0)
               ]
    model = generate_box_model(5, start=10, skip_fibids=[11], missing_relids=[2])

    assert len(model) == len(expected)

    for m, e in zip(model, expected):
        assert m == e


def test_count_peaks1():
    with pytest.raises(ValueError):
        count_peaks([])


def test_count_peaks():

    expected = []
    idx = 0
    for p in PEAKS:
        t = (idx + 1, p, 0, idx)
        expected.append(t)
        idx += 1

    result = count_peaks(PEAKS, tol=1.2, distance=6.0)

    assert result == expected