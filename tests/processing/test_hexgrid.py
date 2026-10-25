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
    hexgrid_outer_points,
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


def test_hexgrid_outer_points_one():
    """The outer points of a point are its 6 neighbors"""
    outer, neighbors = hexgrid_outer_points(np.array([[0.0], [0.0]]))
    assert outer.shape == (2, 6)
    assert np.allclose(np.hypot(outer[0], outer[1]), 1.0)
    assert neighbors == [[0]] * 6


def test_hexgrid_outer_points_lcb():
    """The outer points of the LCB are a ring around it, without repeated points"""
    r0l = np.array(calc_lcb_grid())
    outer, neighbors = hexgrid_outer_points(r0l)
    # the ring of a grid of 27 x 21 points
    assert outer.shape == (2, 2 * 27 + 2 * 21 + 2)
    dist = np.hypot(outer[0][:, np.newaxis] - outer[0], outer[1][:, np.newaxis] - outer[1])
    np.fill_diagonal(dist, np.inf)
    assert dist.min() > 0.99
    for point, nb in zip(outer.T, neighbors):
        assert 1 <= len(nb) <= 3
        assert np.allclose(np.hypot(*(r0l[:, nb] - point[:, np.newaxis])), 1.0)


def test_hexgrid_outer_points_hole():
    """A missing point inside the grid is an outer point"""
    r0l = calc_matrix(5, 5)
    outer, neighbors = hexgrid_outer_points(np.delete(r0l, 12, axis=1))
    inside = [idx for idx, nb in enumerate(neighbors) if len(nb) == 6]
    assert len(inside) == 1
    assert np.allclose(outer[:, inside[0]], r0l[:, 12])


def test_neigh6_level1():
    """The first level are the cell and its 6 connected cells"""
    assert sorted(neigh6(3, 4)) == sorted([(3, 4)] + list(connected6(3, 4)))


@pytest.mark.parametrize("level", [1, 2, 3])
def test_neigh6_size(level):
    """A hexagon of radius level, without repeated cells"""
    cells = list(neigh6(0, 0, level=level))
    assert len(cells) == len(set(cells)) == 3 * level * (level + 1) + 1
