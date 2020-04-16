import pytest
import numpy as np
import math

from ..hexspline import rescaling_kernel, hexspline2

@pytest.mark.parametrize("p", [1, 2])
@pytest.mark.parametrize("scale", [0.5, 0.8, 1])
def test_rescaling_kernel_normalization(p, scale):
    rbs = rescaling_kernel(p, scale=scale)
    expected_area = math.sqrt(3) / 2
    area = rbs.integral(-3.0, 3.0, -3.0, 3.0)
    assert np.allclose(area, expected_area, rtol=1e-2)


def test_hexspline2():
    x = [0.45, -0.2, 1, 0, 0.68, 1.2]
    y = [-0.8, -0.3, 0.5, 0, 0.2, 0.5]
    expected_res = [5.60769515e-02, 5.78341925e-01, 6.40987562e-16, 1.00000000e+00,
       2.11842907e-01, 7.69185075e-16]
    res = hexspline2(x, y)
    assert np.allclose(res, expected_res, rtol=1e-2)