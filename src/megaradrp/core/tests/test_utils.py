
import numpy
import pytest

from ..utils import atleast_2d_last

@pytest.mark.parametrize("arr, shape, ndim", [
    ([1, 2, 3], (3, 1), 2),
    ([1, 2, 3, 789], (4, 1), 2),
    (4, (1, 1), 2),
    ([[1,2,3], [4,5,6]], (2,3), 2),
    (numpy.empty((3,3,3)), (3,3,3), 3)
])
def test_utils_atleast1(arr, shape, ndim):
    # With one input, return de object
    res = atleast_2d_last(arr)
    print(res, res.shape)
    assert res.shape == shape
    assert res.ndim == ndim

@pytest.mark.parametrize("arrs, shapes, ndims", [
    ([[1, 2, 3], [1, 2, 3, 789], 4, [[1, 2, 3], [4, 5, 6]], numpy.empty((3, 3, 3))],
     [(3, 1), (4, 1), (1, 1), (2, 3), (3, 3, 3)],
     [2, 2, 2, 2, 3])
])
def test_utils_atleast2(arrs, shapes, ndims):
    # With more than one, a list
    res = atleast_2d_last(*arrs)
    assert isinstance(res, list)
    for el, shape, ndim in zip(res, shapes, ndims):
        assert el.shape == shape
        assert el.ndim == ndim
