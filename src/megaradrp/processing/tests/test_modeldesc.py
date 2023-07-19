
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

    # import matplotlib.pyplot as plt

    # plt.plot(xval, box)
    # plt.show()

    return x_0, box


def test_0():
    gamma = 3.9
    alpha = 3
    x_0, column = create_column(gamma, alpha)
    model_desc = MoffatModelDescription(np.array(x_0), fixed_center=True)
    values = model_desc.init_values(column)
    assert sorted(values.keys()) == sorted(model_desc.model_cls.param_names)
    for p in model_desc.model_cls.param_names:
        assert len(values[p]) == len(x_0)

