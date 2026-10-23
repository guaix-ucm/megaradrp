import numpy as np
import pytest

from megaradrp.processing.hexgrid import (
    H_HEX,
    M_SQRT3,
    R_HEX,
    calc_lcb_grid,
    calc_matrix,
    connected6,
    hexgrid_extremes,
    neigh6,
)


def test_calc_matrix():
    res = calc_matrix(4, 4, grid_type=2)
    r1 = M_SQRT3 * np.array(
        [
            0,
            0,
            0,
            0,
            1.0 / 2,
            1.0 / 2,
            1.0 / 2,
            1.0 / 2,
            1,
            1,
            1,
            1,
            3.0 / 2,
            3.0 / 2,
            3.0 / 2,
            3.0 / 2,
        ]
    )
    r2 = [0, 1, 2, 3, -0.5, 0.5, 1.5, 2.5, 0, 1, 2, 3, -0.5, 0.5, 1.5, 2.5]
    expected = np.array([r1, r2])
    assert np.allclose(res, expected)


def test_calc_lcb_grid():
    res_x, res_y = calc_lcb_grid()
    assert res_x.shape == (567,)
    assert res_y.shape == (567,)
    assert np.allclose([res_x.min(), res_x.max()], [-11.258330249197702, 11.258330249197702])
    assert np.allclose([res_y.min(), res_y.max()], [-10.5, 10.0])


def test_hexgrid_extremes():
    r0l = np.array([[0.0, 2.0], [0.0, 1.0]])
    assert hexgrid_extremes(r0l, 0.5) == ((-1, 3), (-1, 5))


@pytest.mark.parametrize("target_scale", [0.5, 0.8, 1.0, 2.0])
def test_hexgrid_extremes_cover(target_scale):
    """The rectangular grid covers the hexagons around the points"""
    r0l = np.array(calc_lcb_grid())
    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, target_scale)
    # the borders of the pixels of the rectangular grid
    assert (j1min - 0.5) * target_scale <= r0l[0].min() - R_HEX
    assert (j1max + 0.5) * target_scale >= r0l[0].max() + R_HEX
    assert (i1min - 0.5) * target_scale <= r0l[1].min() - H_HEX
    assert (i1max + 0.5) * target_scale >= r0l[1].max() + H_HEX


def test_neigh6_level1():
    """The first level are the cell and its 6 connected cells"""
    assert sorted(neigh6(3, 4)) == sorted([(3, 4)] + list(connected6(3, 4)))


@pytest.mark.parametrize("level", [1, 2, 3])
def test_neigh6_size(level):
    """A hexagon of radius level, without repeated cells"""
    cells = list(neigh6(0, 0, level=level))
    assert len(cells) == len(set(cells)) == 3 * level * (level + 1) + 1
