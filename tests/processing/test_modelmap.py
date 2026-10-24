from astropy.modeling.functional_models import Moffat1D
import math
import numpy as np
import numpy.polynomial.polynomial as nppol
import pytest
from megaradrp.processing.modelmap import aper_extract, calc1d_model, calc_matrix, calc_matrix_cols
from megaradrp.products.modelmap import GeometricModel, ModelMap
from megaradrp.processing.modeldesc.moffat import MoffatModelDescription
from megaradrp.processing.modeldesc.gaussbox import GaussBoxModelDescription


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
        val = Moffat1D.evaluate(x, x_0=k, amplitude=1, gamma=gamma, alpha=alpha)
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
    valid = range(1, 623 + 1)
    model_desc = MoffatModelDescription(fixed_center=True)
    res = calc1d_model(model_desc, box, x_0, valid, 0, lateral=2, nloop=1)

    ref_values = {
        "x_0": x_0,
        "amplitude": np.ones_like(x_0),
        "gamma": gamma * np.ones_like(x_0),
        "alpha": alpha * np.ones_like(x_0),
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
    # the fibers are numbered from 1
    valid = range(1, 624)
    wshape = (4112, 623)
    params = {}
    params["mean"] = g_mean
    params["stddev"] = g_std
    wm = calc_matrix(wshape, model, params, valid)
    assert isinstance(wm, scipy.sparse.csr_matrix)
    # each fiber is in its column, around its row
    assert wm.shape == wshape
    assert wm[100, 0] > 0.3
    assert wm[106, 1] > 0.3
    assert wm[100, 1] == 0


GAUSSBOX = GaussBoxModelDescription().model_cls


def gaussbox_column(centers, amplitudes, stddev, nrows=4112):
    xl = np.arange(nrows, dtype=float)
    column = np.zeros(nrows)
    for center, amplitude in zip(centers, amplitudes):
        column += GAUSSBOX.evaluate(xl, amplitude, center, stddev, 0.5)
    return column


def test_calc1d_model_gaussbox():
    """The fit recovers the width and amplitude of the profiles"""
    np.random.seed(1)
    valid = list(range(101, 111))
    centers = 6.33 * np.array(valid, dtype=float)
    amplitudes = np.linspace(80.0, 120.0, len(valid))
    column = gaussbox_column(centers, amplitudes, 1.7)

    model_desc = GaussBoxModelDescription(sigma=1.53)
    result = calc1d_model(model_desc, column, centers, valid, 2000, lateral=2, nloop=2)

    assert sorted(result) == valid
    stddev = np.array([result[fibid]["stddev"] for fibid in valid])
    amplitude = np.array([result[fibid]["amplitude"] for fibid in valid])
    assert np.allclose(stddev, 1.7, rtol=1e-3)
    assert np.allclose(amplitude, amplitudes, rtol=1e-3)
    # the centers are fixed
    assert np.allclose([result[fibid]["mean"] for fibid in valid], centers)


def test_calc1d_model_shape():
    column = np.zeros(100)
    with pytest.raises(ValueError):
        calc1d_model(GaussBoxModelDescription(), column, [10.0], [1], 0)


def make_model_map(means, stddev, valid=None, offset=0.0):
    """A ModelMap of gaussbox profiles, constant along the columns"""
    model_map = ModelMap(instrument="MEGARA")
    model_map.total_fibers = len(means)
    model_map.global_offset = nppol.Polynomial([offset])
    model_map.ref_column = 0
    for fibid, mean in enumerate(means, 1):
        if valid is None or fibid in valid:
            params = {"mean": nppol.Polynomial([mean]), "stddev": nppol.Polynomial([stddev])}
            model = {"model_name": "gaussbox", "params": params}
        else:
            model = {}
        model_map.contents.append(GeometricModel(fibid, 1, 0, 7, model))
    return model_map


def make_image(means, stddev, fluxes, nrows=120):
    """Image with the profiles of the fibers, with a flux per fiber and column"""
    rows = np.arange(nrows, dtype=float)
    img = np.zeros((nrows, fluxes.shape[1]))
    for mean, flux in zip(means, fluxes):
        profile = GAUSSBOX.evaluate(rows, 1.0, mean, stddev, 0.5)
        img += profile[:, np.newaxis] * flux
    return img


@pytest.mark.parametrize("offset", [0.0, 1.5])
def test_aper_extract(offset):
    """The extraction with the model recovers the flux of overlapping fibers"""
    means = 20.0 + 6.0 * np.arange(6)
    fluxes = np.outer(np.arange(1.0, 7.0), np.linspace(1.0, 2.0, 8)) * 100
    img = make_image(means + offset, 1.6, fluxes)
    model_map = make_model_map(means, 1.6, offset=offset)

    wcols = calc_matrix_cols(model_map, img.shape)
    assert sorted(wcols) == list(range(8))
    rss = aper_extract(model_map, wcols, img)
    assert rss.shape == (6, 8)
    assert np.allclose(rss, fluxes, rtol=1e-4)


def test_aper_extract_invalid_fiber():
    """The invalid fibers are not extracted"""
    means = 20.0 + 6.0 * np.arange(6)
    fluxes = np.full((6, 8), 100.0)
    fluxes[2] = 0.0
    img = make_image(means, 1.6, fluxes)
    model_map = make_model_map(means, 1.6, valid=[1, 2, 4, 5, 6])

    rss = aper_extract(model_map, calc_matrix_cols(model_map, img.shape), img)
    assert np.all(rss[2] == 0)
    assert np.allclose(rss[[0, 1, 3, 4, 5]], 100.0, rtol=1e-4)


def test_calc_matrix_cols_undefined_model():
    model_map = make_model_map([20.0], 1.6)
    model_map.contents[0].model["model_name"] = "undefined"
    with pytest.raises(ValueError):
        calc_matrix_cols(model_map, (120, 8))


@pytest.mark.slow
@pytest.mark.xfail(reason="the pool of processes is not closed, ResourceWarning", strict=True)
def test_calc_matrix_cols_parallel():
    means = 20.0 + 6.0 * np.arange(6)
    fluxes = np.full((6, 8), 100.0)
    img = make_image(means, 1.6, fluxes)
    model_map = make_model_map(means, 1.6)
    wcols = calc_matrix_cols(model_map, img.shape, processes=2)
    assert np.allclose(aper_extract(model_map, wcols, img), 100.0, rtol=1e-4)
