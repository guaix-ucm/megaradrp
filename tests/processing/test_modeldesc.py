from astropy.modeling.functional_models import Moffat1D
import math
import numpy as np
import pytest
from megaradrp.processing.modeldesc.moffat import MoffatModelDescription
from megaradrp.processing.modeldesc.gaussbox import GaussBoxModelDescription


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
        eval = Moffat1D.evaluate(x, x_0=k, amplitude=1, gamma=gamma, alpha=alpha)
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

    kwds = {"npix": 8}
    model_desc = model_cls.create_from_centers(centers, **kwds)

    assert isinstance(model_desc, GaussBoxModelDescription)


def gaussian_column(centers, amplitude, sigma, size=200):
    xval = np.arange(size, dtype=float)
    column = np.zeros(size)
    for center in centers:
        column += amplitude * np.exp(-((xval - center) ** 2) / (2 * sigma**2))
    return column


def test_gaussbox_create_from_centers():
    model_desc = GaussBoxModelDescription.create_from_centers(fixed_center=False, npix=7, sigma=2.0, init_simple=True)
    assert model_desc.fixed_center is False
    assert (model_desc.npix, model_desc.sigma, model_desc.init_simple) == (7, 2.0, True)
    assert model_desc.params_save == ["mean", "stddev"]
    assert model_desc.deg_save == [3, 5]


def test_gaussbox_init_values():
    centers = np.array([40.0, 80.3, 120.7])
    column = gaussian_column(centers, 10.0, 1.5)
    values = GaussBoxModelDescription(sigma=1.5).init_values(column, centers)
    assert np.allclose(values["mean"], centers)
    assert np.allclose(values["stddev"], 1.5)
    assert np.allclose(values["hpix"], 0.5)
    # the value at the pixel of the center, divided by 0.25
    assert np.allclose(values["amplitude"], column[[40, 80, 121]] / 0.25)


def test_gaussbox_init_values_simple():
    """The width of Gaussian profiles is fitted with a parabola in the logarithm"""
    centers = np.array([40.0, 80.0, 120.0])
    column = gaussian_column(centers, 10.0, 1.5)
    values = GaussBoxModelDescription(init_simple=True, npix=5).init_values(column, centers)
    assert np.allclose(values["mean"], centers)
    assert np.allclose(values["stddev"], 1.5)


@pytest.mark.xfail(reason="the amplitude of init_simple uses mu = b / sig2 and exp(a - mu**2)", strict=True)
def test_gaussbox_init_values_simple_amplitude():
    centers = np.array([40.3, 80.3, 120.3])
    column = gaussian_column(centers, 10.0, 1.5)
    values = GaussBoxModelDescription(init_simple=True, npix=5).init_values(column, centers)
    assert np.allclose(values["amplitude"], 10.0)


@pytest.mark.parametrize("fixed_center", [True, False])
def test_gaussbox_params_fixed_bounds(fixed_center):
    model_desc = GaussBoxModelDescription(fixed_center=fixed_center)
    values = {"mean": 50.0, "stddev": 1.5, "amplitude": 3.0, "hpix": 0.5}

    fixed = model_desc.params_fixed(values, 1, 100)
    assert fixed["hpix"] is True
    assert fixed.get("mean", False) is fixed_center

    bounds = model_desc.params_bounds(values, 1, 100)
    assert bounds["stddev"] == (1.0, 2.0)
    assert bounds["amplitude"] == (0, None)
    if fixed_center:
        assert "mean" not in bounds
    else:
        assert bounds["mean"] == (49.5, 50.5)

    assert model_desc.fiber_center(values) == 50.0


def test_gaussbox_init_values_per_profile():
    centers = np.array([40.0, 80.0])
    column = gaussian_column(centers, 10.0, 1.5)
    model_desc = GaussBoxModelDescription(sigma=1.5)
    per_profile = model_desc.init_values_per_profile(column, centers, [7, 9])
    assert sorted(per_profile) == [7, 9]
    assert per_profile[9]["mean"] == 80.0
    assert sorted(per_profile[7]) == sorted(model_desc.model_cls.param_names)
