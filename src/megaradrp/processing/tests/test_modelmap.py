
from astropy.modeling.functional_models import Moffat1D
import math
import numpy as np
import pytest
from ..modelmap import calc1d_M
from ..modeldesc import MoffatModelDescription


def create_column(gamma, alpha):
    xval = np.arange(4112, dtype=float)
    box = np.zeros((4112,), dtype=float)
    x_0 = []
    for i in range(623):
        k = 6.33 * (i + 1)
        k0 = int(math.floor(k - 5))
        k1 = int(math.ceil(k + 5))
        x_0.append(k)
        x = xval[k0:k1]
        eval = Moffat1D.evaluate(x, x_0=k, amplitude=1, gamma=gamma, alpha=alpha)
        box[k0:k1] += eval
    return x_0, box


def test_1():
    gamma = 3.9
    alpha = 3
    x_0, y_col = create_column(gamma, alpha)
    valid = range(1, 623+1)
    model_desc = MoffatModelDescription(np.array(x_0), fixed_center=True)
    res = calc1d_M(model_desc, y_col, valid, 0, lateral=2, nloop=1)

    ref_values = {
        'x_0': x_0,
        'amplitude': np.ones_like(x_0),
        'gamma': gamma * np.ones_like(x_0),
        'alpha': alpha * np.ones_like(x_0),
    }

    for fibid in valid:
        fiber_fit = res.get(fibid)
        # The fit is done
        assert fiber_fit is not None

        for param in Moffat1D.param_names:
            assert param in fiber_fit
            print(fibid, param, fiber_fit[param], ref_values[param][fibid - 1])


def test_2():
    box = np.zeros((4112,), dtype=float)
    valid = [2, 12, 34]
    helper = MoffatModelDescription(np.array([3, 38]))
    with pytest.raises(ValueError):
        calc1d_M(helper, box, valid, 0)
