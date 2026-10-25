import pytest
import numpy as np
import math

from scipy import signal
import scipy.interpolate as interp

from megaradrp.processing.hexspline import HEXSPLINE3_LATTICE, hexspline1, hexspline2, rescaling_kernel
from megaradrp.simulation.convolution import setup_grid


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


def test_hexspline3_lattice():
    """The values in the grid of the hexspline of order 3, tabulated by convolution"""
    step = 0.005
    xx, yy, xs, ys, _, _ = setup_grid(4.0, 4.0, step, step)
    omega = math.sqrt(3) / 2
    hex1 = hexspline1(xx, yy)
    hex2 = signal.fftconvolve(hex1, hex1, mode="same") * step**2 / omega
    hex3 = signal.fftconvolve(hex2, hex1, mode="same") * step**2 / omega
    rbs = interp.RectBivariateSpline(xs, ys, hex3)
    # (x, y) of the 6 neighbors
    angles = np.radians([30, 90, 150, 210, 270, 330])
    eta0, eta1 = HEXSPLINE3_LATTICE
    assert math.isclose(eta0 + 6 * eta1, 1.0)
    assert np.allclose(rbs.ev(0.0, 0.0), eta0, atol=1e-4)
    assert np.allclose(rbs.ev(np.sin(angles), np.cos(angles)).mean(), eta1, atol=1e-4)
    # zero in the second neighbors, in the border of the support
    angles2 = np.radians([0, 60, 120, 180, 240, 300])
    assert np.allclose(rbs.ev(math.sqrt(3) * np.sin(angles2), math.sqrt(3) * np.cos(angles2)), 0.0, atol=1e-4)
