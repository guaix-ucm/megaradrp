import astropy.io.fits as fits
import numpy
import pytest
from numina.core import DataFrame
from numina.modeling.gaussbox import GaussBox

import megaradrp.recipes.calibration.modelmap as modelmap
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.products.tracemap import GeometricTrace, TraceMap
from megaradrp.recipes.calibration.modelmap import ModelMapRecipe, column_to_fit
from megaradrp.testing.common import make_obresult, raw_header

SIGMA = 1.53
# fibers with a trace, the others are invalid
VALID = range(101, 111)


def center(fibid, x):
    """The center of the profile of the fiber, in the 0-based column x"""
    return 1000.0 + 7.0 * (fibid - 101) + 2e-4 * x


def amplitude(x):
    return 1000.0 * (1.0 + 0.2 * numpy.sin(x / 1000.0))


@pytest.fixture(scope="module")
def reduced():
    """A reduced image with the profiles of the valid fibers"""
    rows = numpy.arange(4112)[:, numpy.newaxis]
    cols = numpy.arange(4096)
    data = numpy.zeros((4112, 4096), dtype="float32")
    for fibid in VALID:
        sl = slice(int(center(fibid, 0)) - 10, int(center(fibid, 4096)) + 10)
        data[sl] += GaussBox.evaluate(rows[sl], amplitude(cols), center(fibid, cols), SIGMA, 0.5)
    hdr = raw_header(vph="LR-R", insmode="LCB")
    hdr["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    fibers = fits.ImageHDU(header=create_default_fiber_header("LCB"), name="FIBERS")
    return fits.HDUList([fits.PrimaryHDU(data, header=hdr), fibers])


def make_tracemap():
    tracemap = TraceMap(instrument="MEGARA")
    tracemap.total_fibers = 623
    tracemap.missing_fibers = [623]
    for fibid in range(1, 624):
        fitparms = [center(fibid, 0), 2e-4] if fibid in VALID else []
        tracemap.contents.append(GeometricTrace(fibid, 1, 1, 4096, fitparms=fitparms))
    return tracemap


@pytest.fixture
def run_recipe(monkeypatch, reduced):
    """Run the recipe, with the basic reduction replaced by the reduced image"""

    def run(recipe, **kwds):
        monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: None)
        monkeypatch.setattr(modelmap, "basic_processing_with_combination", lambda rinput, flow, method: reduced)
        obresult = make_obresult("MegaraModelMap", header=raw_header(vph="LR-R", insmode="LCB"))
        raw = fits.PrimaryHDU(header=raw_header(vph="LR-R", insmode="LCB"))
        raw.header["SENTEMP4"] = 290.15
        raw.header["UUID"] = "11111111-2222-3333-4444-555555555555"
        obresult.frames = [DataFrame(frame=fits.HDUList([raw]))]
        # the bias is not used, the reduction is replaced
        master_bias = DataFrame(frame=fits.HDUList([fits.PrimaryHDU()]))
        values = dict(processes=1, debug_plot=0, master_bias=master_bias)
        values.update(kwds)
        rinput = recipe.create_input(obresult=obresult, master_traces=make_tracemap(), **values)
        return recipe.run(rinput)

    return run


@pytest.mark.slow
def test_modelmap(run_recipe):
    result = run_recipe(ModelMapRecipe())
    model_map = result.master_model
    assert model_map.total_fibers == 623
    assert model_map.missing_fibers == [623]
    assert model_map.tags["vph"] == "LR-R"
    assert model_map.tags["insmode"] == "LCB"
    assert model_map.tags["temp"] == 17.0
    assert [model.fibid for model in model_map.contents] == list(VALID)

    cols = numpy.arange(200, 3900, 100)
    for model in model_map.contents:
        assert model.model["model_name"] == "gaussbox"
        params = model.model["params"]
        assert numpy.allclose(params["mean"](cols), center(model.fibid, cols), atol=0.02)
        assert numpy.allclose(params["stddev"](cols), SIGMA, atol=0.02)

    # the extraction with the model recovers the flux of the fibers
    rss = result.reduced_rss.open()[0].data
    rows = numpy.arange(4112)
    for fibid in VALID:
        flux = GaussBox.evaluate(rows, amplitude(2000), center(fibid, 2000), SIGMA, 0.5).sum()
        assert rss[fibid - 1, 2000] == pytest.approx(flux, rel=0.002)
    assert numpy.all(rss[0] == 0)


@pytest.mark.slow
def test_modelmap_intermediate(run_recipe, tmp_path, monkeypatch):
    """The plots of the parameters and the DS9 regions are saved"""
    monkeypatch.chdir(tmp_path)
    run_recipe(ModelMapRecipe(intermediate_results=True), debug_plot=1)
    names = {path.name for path in tmp_path.iterdir()}
    assert {"reduced_image.fits", "ds9.reg", "ds9_raw.reg"} <= names
    assert {f"fib_{fibid:03d}_{name}.png" for fibid in VALID for name in ["mean", "stddev"]} <= names


def test_modelmap_processes(run_recipe, monkeypatch):
    """With processes=0, the number of processes depends on the CPUs"""
    used = []

    def fit_model(model_obj, data, tracemap, cols, processes):
        used.append(processes)
        raise RuntimeError("stop")

    monkeypatch.setattr(modelmap, "fit_model", fit_model)
    monkeypatch.setattr(modelmap.mp, "cpu_count", lambda: 8)
    with pytest.raises(RuntimeError):
        run_recipe(ModelMapRecipe(), processes=0)
    monkeypatch.setattr(modelmap.mp, "cpu_count", lambda: 2)
    with pytest.raises(RuntimeError):
        run_recipe(ModelMapRecipe(), processes=0)
    assert used == [6, 1]


def test_column_to_fit():
    data = numpy.arange(10.0) * numpy.ones((3, 1))
    assert numpy.all(column_to_fit(data, 5) == 5)


def test_column_to_fit_average():
    """The column is averaged with the columns around it"""
    data = numpy.arange(10.0) ** 2 * numpy.ones((3, 1))
    expected = numpy.mean([3**2, 4**2, 5**2, 6**2, 7**2])
    assert numpy.allclose(column_to_fit(data, 5, average=2), expected)
