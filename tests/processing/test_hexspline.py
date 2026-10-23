import pytest
import numpy as np
import math

from megaradrp.processing.hexspline import rescaling_kernel, hexspline2


# scales of the cubes with pixels of 0.3, 0.5 and 1 arcsec
@pytest.mark.parametrize(
    "p, scale",
    [
        (1, 0.56),
        (1, 0.93),
        (1, 1.86),
        (2, 0.56),
        (2, 0.93),
        pytest.param(2, 1.86, marks=pytest.mark.slow),
        pytest.param(3, 1.0, marks=pytest.mark.slow),
    ],
)
def test_rescaling_kernel_normalization(p, scale):
    rbs = rescaling_kernel(p, scale=scale)
    expected_area = math.sqrt(3) / 2
    area = rbs.integral(-10.0, 10.0, -10.0, 10.0)
    assert np.allclose(area, expected_area, rtol=1e-2)


@pytest.mark.parametrize("p, scale", [(1, 0.56), (1, 1.86), (2, 0.56), pytest.param(2, 1.86, marks=pytest.mark.slow)])
def test_rescaling_kernel_support(p, scale):
    """The kernel is zero far from the origin (#343)"""
    rbs = rescaling_kernel(p, scale=scale)
    far = [5.0, 10.0, 20.0]
    assert np.allclose(rbs.ev(far, [0.0, 0.0, 0.0]), 0.0, atol=1e-8)
    assert np.allclose(rbs.ev([0.0, 0.0, 0.0], far), 0.0, atol=1e-8)


def test_hexspline2():
    x = [0.45, -0.2, 1, 0, 0.68, 1.2]
    y = [-0.8, -0.3, 0.5, 0, 0.2, 0.5]
    expected_res = [
        5.60769515e-02,
        5.78341925e-01,
        6.40987562e-16,
        1.00000000e00,
        2.11842907e-01,
        7.69185075e-16,
    ]
    res = hexspline2(x, y)
    assert np.allclose(res, expected_res, rtol=1e-2)


def test_rescaling_kernel_raise():
    with pytest.raises(ValueError):
        rescaling_kernel(4)
