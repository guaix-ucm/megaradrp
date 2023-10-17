
from astropy.modeling.functional_models import Moffat1D
import math
import numpy as np
from ..modeldesc.moffat import MoffatModelDescription
from ..modeldesc.gaussbox import GaussBoxModelDescription


def create_moffat_column(gamma, alpha):
    xval = np.arange(4112, dtype=float)
    box = np.zeros((4112,), dtype=float)
    x_0 = []
    for i in range(623):
        k = 6.33 * (i + 1)
        k0 = int(math.floor(k - 5))
        k1 = int(math.ceil(k + 5))
        x_0.append(k)
        x = xval[k0:k1]
        eval = Moffat1D.evaluate(
            x, x_0=k, amplitude=1, gamma=gamma, alpha=alpha)
        box[k0:k1] += eval

    return x_0, box


def test_0():
    gamma = 3.9
    alpha = 3
    x_0, column = create_moffat_column(gamma, alpha)
    model_desc = MoffatModelDescription(fixed_center=True)
    values = model_desc.init_values(column, x_0)
    assert sorted(values.keys()) == sorted(model_desc.model_cls.param_names)
    for p in model_desc.model_cls.param_names:
        assert len(values[p]) == len(x_0)


def test_3():
    centers = np.array([4, 5, 6])

    model_cls = GaussBoxModelDescription

    kwds = {'npix': 8}
    model_desc = model_cls.create_from_centers(centers, **kwds)

    assert isinstance(model_desc, GaussBoxModelDescription)
