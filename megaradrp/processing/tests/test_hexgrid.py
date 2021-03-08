
import numpy as np
from ..hexgrid import calc_matrix, M_SQRT3, calc_lcb_grid


def test_calc_matrix():
    res = calc_matrix(4, 4, grid_type=2)
    r1 = M_SQRT3 * np.array(
        [0, 0, 0, 0, 1.0/2, 1.0/2, 1.0/2, 1.0/2, 1, 1, 1, 1, 3.0/2, 3.0/2, 3.0/2, 3.0/2]
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
