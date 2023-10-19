
from astropy.modeling.functional_models import Moffat1D
import math
import numpy as np
import pytest
from ..modelmap import calc1d_model, calc_matrix
from ..modeldesc.moffat import MoffatModelDescription
from ..modeldesc.gaussbox import GaussBoxModelDescription


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
        val = Moffat1D.evaluate(
            x, x_0=k, amplitude=1, gamma=gamma, alpha=alpha
        )
        box[k0:k1] += val

    # import matplotlib.pyplot as plt

    # plt.plot(xval, box)
    # plt.show()

    return x_0, box


def test_0():
    gamma = 3.9
    alpha = 3
    x_0, column = create_column(gamma, alpha)
    # valid = range(1, 623 + 1)
    model_desc = MoffatModelDescription(fixed_center=True)
    values = model_desc.init_values(column, x_0)
    assert sorted(values.keys()) == sorted(model_desc.model_cls.param_names)


def _test_1():
    gamma = 3.9
    alpha = 3
    x_0, box = create_column(gamma, alpha)
    valid = range(1, 623+1)
    model_desc = MoffatModelDescription(fixed_center=True)
    res = calc1d_model(model_desc, box, x_0, valid, 0, lateral=2, nloop=1)

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
    centers = np.array([3, 38])
    helper = MoffatModelDescription()
    with pytest.raises(ValueError):
        calc1d_model(helper, box, centers, valid, 0)


def test_calc1():
    import scipy.sparse
    model = GaussBoxModelDescription().model_cls
    g_mean = 100 + 6 * np.arange(623)
    g_std = 1 + np.zeros_like(g_mean)
    valid = range(623)
    wshape = (4112, 623)
    params = {}
    params['mean'] = g_mean
    params['stddev'] = g_std
    wm = calc_matrix(wshape, model, params, valid)
    assert isinstance(wm, scipy.sparse.csr_matrix)
