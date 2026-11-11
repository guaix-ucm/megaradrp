import types

from astropy.io import fits
import numpy
import pytest
import numina.exceptions
from numina.types.qc import QC
from numina.util.node import Node

import megaradrp.recipes.calibration.flat as flat
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.recipes.calibration.flat import (
    FiberFlatRecipe,
    _flat_knots,
    _trace_polynomials,
    _smoothing_window_check,
    absorption_windows,
    check_qc_fiberflat,
    find_zero_order,
    flag_fibers,
    interpolate_windows,
    remove_zero_order,
)
from megaradrp.testing.create_image import create_rss

NWL = 4300


@pytest.fixture(scope="module")
def transmission():
    return numpy.random.default_rng(0).uniform(0.8, 1.2, 623)


def flat_rss(transmission, lamp, invalid=()):
    """Calibrated RSS of a flat: the spectrum of the lamp times the transmission of each fiber"""
    data = (transmission[:, numpy.newaxis] * lamp).astype("float32")
    rss = create_rss(data, numpy.ones_like(data, dtype="int16"))
    hdrf = create_default_fiber_header("LCB")
    for fibid in invalid:
        hdrf[f"FIB{fibid:03d}_V"] = False
    rss["FIBERS"].header = hdrf
    return rss


def expected_flat(transmission, invalid=()):
    valid = [idx for idx in range(623) if idx + 1 not in invalid]
    return valid, transmission / transmission[valid].mean()


def test_obtain_fiber_flat(transmission):
    """With a flat lamp, the fiber flat is the relative transmission of each fiber"""
    rss = flat_rss(transmission, numpy.ones(NWL))
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=1e-4)
    # the extensions are copied
    assert "WLMAP" in result
    assert result is not rss


def test_obtain_fiber_flat_no_coverage(transmission):
    """The pixels without coverage are 1"""
    rss = flat_rss(transmission, numpy.ones(NWL))
    rss["WLMAP"].data[5, :100] = 0
    rss[0].data[5, :100] = 0
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    assert numpy.all(result[0].data[5, :100] == 1.0)


def test_obtain_fiber_flat_invalid_fibers(transmission):
    """The invalid fibers do not change the normalization"""
    invalid = list(range(1, 31))
    rss = flat_rss(transmission, numpy.ones(NWL), invalid=invalid)
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission, invalid)
    assert numpy.allclose(result[0].data[valid], expected[valid, numpy.newaxis], rtol=1e-3)


def test_obtain_fiber_flat_lamp_shape(transmission):
    """The spectrum of the lamp is removed"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    rss = flat_rss(transmission, lamp)
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=1e-2)


@pytest.mark.parametrize("value, valid", [(31, True), (5, True), (3, False), (32, False)])
def test_smoothing_window_check(value, valid):
    if valid:
        assert _smoothing_window_check(value) == value
    else:
        with pytest.raises(numina.exceptions.ValidationError):
            _smoothing_window_check(value)


def test_process_flat2d(monkeypatch):
    recipe = FiberFlatRecipe()
    calls = []
    img = fits.HDUList([fits.PrimaryHDU(numpy.ones((10, 10)))])

    def combination(rinput, flow, method, method_kwargs):
        calls.append((flow, method.__name__, method_kwargs))
        return img

    monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: "flow")
    monkeypatch.setattr(flat, "basic_processing_with_combination", combination)
    rinput = types.SimpleNamespace(
        obresult=types.SimpleNamespace(configuration=None), method="median", method_kwargs={"a": 1}
    )
    result = recipe.process_flat2d(rinput)
    assert calls == [("flow", "median", {"a": 1})]
    assert result is img
    assert "NUMXVER" in img[0].header


def test_fiber_flat_recipe(monkeypatch, transmission):
    """The recipe returns the 2D image and the RSS before the flip and the wavelength calibration"""
    calls = []
    img2d = fits.HDUList([fits.PrimaryHDU(numpy.full((10, 10), 7.0))])
    raw_rss = flat_rss(transmission, numpy.ones(NWL))
    raw_rss[0].data[:, :10] = 3.0

    class Extractor(Node):
        def __init__(self, tracemap, datamodel, offset):
            super().__init__()
            calls.append(("extractor", tracemap, offset))

        def run(self, img):
            return raw_rss

    class Calibrator(Node):
        def __init__(self, wlcalib, datamodel):
            super().__init__()
            calls.append(("calibrator", wlcalib))

        def run(self, img):
            calls.append(("calibrate", float(img[0].data[0, -1])))
            return img

    recipe = FiberFlatRecipe()
    saved = []
    monkeypatch.setattr(recipe, "process_flat2d", lambda rinput: img2d)
    monkeypatch.setattr(recipe, "save_intermediate_img", lambda img, name: saved.append(name))
    monkeypatch.setattr(flat, "ApertureExtractor", Extractor)
    monkeypatch.setattr(flat, "WavelengthCalibrator", Calibrator)

    rinput = types.SimpleNamespace(
        master_apertures="traces",
        master_wlcalib="wlcalib",
        extraction_offset=[0.5],
        smoothing_window=31,
        smoothing_knot_step=150.0,
        skip_absorption_windows=True,
        remove_zero_order=False,
    )
    result = recipe.run(rinput)

    assert saved == ["reduced_image.fits", "reduced_rss.fits"]
    # the wavelength calibration receives the flipped RSS
    assert calls == [("extractor", "traces", [0.5]), ("calibrator", "wlcalib"), ("calibrate", 3.0)]
    assert numpy.allclose(result.reduced_image.open()[0].data, 7.0)
    # the reduced RSS is the copy before the flip
    assert numpy.allclose(result.reduced_rss.open()[0].data[:, :10], 3.0)
    master = result.master_fiberflat.open()
    assert master[0].header["NUMTYPE"] == "MasterFiberFlat"
    assert master[0].data.shape == (623, NWL)


def test_obtain_fiber_flat_intermediate(transmission, tmp_path, monkeypatch):
    """The plot of the collapsed spectrum is saved without pyplot"""
    import matplotlib.pyplot as plt

    monkeypatch.chdir(tmp_path)
    rss = flat_rss(transmission, numpy.ones(NWL))
    before = plt.get_fignums()
    FiberFlatRecipe(intermediate_results=True).obtain_fiber_flat(rss)
    names = {path.name for path in tmp_path.iterdir()}
    assert names == {"collapse.txt", "mask_noinfo.txt", "collapsed_smooth.png"}
    assert plt.get_fignums() == before


def test_obtain_fiber_flat_faint_fiber(transmission):
    """A fiber much fainter than the others, with a wrong spectrum, does not change the flat of the others"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    rss = flat_rss(transmission, lamp)
    rss[0].data[100] = 0.05 * (1 + (x > 3000) * 50.0)
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission, invalid=[101])
    others = [idx for idx in valid if idx != 100]
    assert numpy.allclose(result[0].data[others], expected[others, numpy.newaxis], rtol=1e-2)


def test_obtain_fiber_flat_partial_coverage(transmission):
    """Near the borders only some fibers have coverage, the average is divided by the fibers used"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    rss = flat_rss(transmission, lamp)
    # the first 300 fibers start 100 pixels later
    rss["WLMAP"].data[:300, :100] = 0
    rss[0].data[:300, :100] = 0
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data[300:, :100], expected[300:, numpy.newaxis], rtol=1e-2)
    assert numpy.all(result[0].data[:300, :100] == 1.0)


def test_obtain_fiber_flat_legacy_spline(transmission):
    """With knot_step <= 0, the spline of degree 5 with the default smoothing is used"""
    rss = flat_rss(transmission, numpy.ones(NWL))
    result = FiberFlatRecipe().obtain_fiber_flat(rss, knot_step=0)
    valid, expected = expected_flat(transmission)
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=1e-4)


def test_absorption_windows():
    """The windows of each resolution class include the O2 A band and the Na D lines"""
    for resolution in ("LR", "MR", "HR"):
        windows = absorption_windows(resolution)
        assert any(w1 <= 7600 <= w2 for w1, w2 in windows)
        assert any(w1 <= 5890 <= w2 for w1, w2 in windows)
        assert all(w1 < w2 for w1, w2 in windows)
    assert absorption_windows("XX") == []


def test_interpolate_windows():
    """The skipped pixels are a straight line between both sides; pixels without coverage are not used"""
    x = numpy.arange(100, dtype=float)
    data = numpy.vstack([1.0 + 0.01 * x, 2.0 + 0.0 * x])
    data[:, 40:50] = 10.0
    skip = numpy.zeros(100, dtype=bool)
    skip[40:50] = True
    coverage = numpy.ones_like(data, dtype=bool)
    interpolate_windows(data, skip, coverage)
    assert numpy.allclose(data[0, 40:50], 1.0 + 0.01 * x[40:50], atol=0.01)
    assert numpy.allclose(data[1, 40:50], 2.0)


def test_obtain_fiber_flat_skip_windows(transmission):
    """An absorption common to all the fibers, inside a window, is not in the fiber flat"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    lamp *= 1 - 0.03 * numpy.exp(-0.5 * ((x - 3000) / 2.0) ** 2)
    rss = flat_rss(transmission, lamp)
    rss[0].header["CRVAL1"] = 7000.0
    rss[0].header["CDELT1"] = 0.4
    rss[0].header["CRPIX1"] = 1.0
    wl_line = 7000.0 + 0.4 * 3000
    valid, expected = expected_flat(transmission)
    result = FiberFlatRecipe().obtain_fiber_flat(rss, skip_windows=[(wl_line - 5, wl_line + 5)])
    assert numpy.allclose(result[0].data, expected[:, numpy.newaxis], rtol=2e-3)
    # without the window, the line stays in the flat
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    assert numpy.min(result[0].data[:, 3000] / expected) < 0.98


def test_flat_knots():
    """The knots are equally spaced, not farther than the step, and need 10 fitted pixels since the previous one"""
    x_fit = numpy.arange(100, 1100)
    knots = _flat_knots(x_fit, 300)
    assert numpy.allclose(knots, [349.75, 599.5, 849.25])
    # a gap without fitted pixels removes the knot after it
    x_fit = numpy.r_[numpy.arange(100, 340), numpy.arange(598, 1100)]
    assert len(_flat_knots(x_fit, 300)) == 2


def test_obtain_fiber_flat_edge_extrapolated(transmission):
    """In the zone of partial coverage the smoothing is continued in a straight line

    The fibers that reach the border are 5 % fainter there than the others: their
    fiber flat keeps it, it is not absorbed by the average of the few fibers there.
    """
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.2 * x / NWL
    rss = flat_rss(transmission, lamp)
    # each group of fibers starts later, from 0 to 150 pixels
    starts = (numpy.arange(623) % 16) * 10
    response = 1.0 - 0.05 * numpy.exp(-x / 100.0)
    for idx, start in enumerate(starts):
        rss["WLMAP"].data[idx, :start] = 0
        rss[0].data[idx, :start] = 0
        if start == 0:
            rss[0].data[idx] *= response
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    valid, expected = expected_flat(transmission)
    ratio = result[0].data[:, :150] / expected[:, numpy.newaxis]
    early = starts == 0
    assert numpy.allclose(ratio[early], response[:150], atol=5e-3)
    covered = rss["WLMAP"].data[~early, :150] > 0
    assert numpy.allclose(ratio[~early][covered], 1.0, atol=5e-3)


def test_obtain_fiber_flat_flags(transmission):
    """Faint and deviant fibers are flagged and stay valid, a fiber without flux is marked invalid"""
    x = numpy.arange(NWL)
    lamp = 1.0 + 0.5 * numpy.exp(-(((x - 2000) / 800.0) ** 2))
    rss = flat_rss(transmission, lamp)
    rss[0].data[100] *= 0.1
    # the response of the fiber 201 falls to the blue and to the red
    rss[0].data[200] *= numpy.exp(-(((x - 2000) / 1200.0) ** 2))
    rss[0].data[300] = 0.0
    result = FiberFlatRecipe().obtain_fiber_flat(rss)
    hdr = result["FIBERS"].header
    assert hdr["FIB101QF"] == "faint"
    assert hdr.get("FIB101_V", True)
    assert hdr["FIB201QF"] == "deviant"
    assert hdr.get("FIB201_V", True)
    assert hdr["FIB301QF"] == "broken"
    assert not hdr["FIB301_V"]
    assert "FIB050QF" not in hdr
    assert check_qc_fiberflat(hdr) == QC.PARTIAL


def test_flag_fibers():
    flat = numpy.ones((12, 400))
    coverage = numpy.ones_like(flat, dtype=bool)
    faint = numpy.zeros(12, dtype=bool)
    faint[1] = True
    good = numpy.ones(12, dtype=bool)
    good[11] = False
    # same shape as the neighbors, other level
    flat[1] = 0.5
    # a step of 50 % in 100 pixels
    flat[5, 300:] = 1.5
    # no flux
    flat[8] = 0.0
    # invalid fibers are not checked
    flat[11] = 7.0
    flat[11, 300:] = 0
    flags = flag_fibers(flat, coverage, faint, good, 10, 20)
    assert flags == {2: ["faint"], 6: ["deviant"], 9: ["broken"]}
    # a step of 20 % is the color of the fiber
    flat[5, 300:] = 1.2
    assert flag_fibers(flat, coverage, faint, good, 10, 20) == {2: ["faint"], 9: ["broken"]}
    # the deviant pixels inside the windows are not counted
    flat[5, 300:] = 1.5
    skip = numpy.zeros(400, dtype=bool)
    skip[290:] = True
    assert flag_fibers(flat, coverage, faint, good, 10, 20, skip=skip) == {2: ["faint"], 9: ["broken"]}
    # non finite values with coverage
    flat[3, 50] = numpy.nan
    assert flag_fibers(flat, coverage, faint, good, 10, 20, skip=skip)[4] == ["broken"]


@pytest.mark.parametrize("ninvalid, qc", [(0, QC.GOOD), (3, QC.PARTIAL), (40, QC.BAD)])
def test_check_qc_fiberflat(ninvalid, qc):
    hdr = create_default_fiber_header("LCB")
    hdr["FIB300QF"] = "faint"
    hdr["FIB301QF"] = "deviant"
    for fibid in range(1, ninvalid + 1):
        hdr[f"FIB{fibid:03d}QF"] = "faint,broken"
    assert check_qc_fiberflat(hdr) == qc


def test_fiber_flat_run_qc(transmission):
    rss = flat_rss(transmission, numpy.ones(NWL))
    rss["FIBERS"].header["FIB010QF"] = "broken"
    recipe = FiberFlatRecipe()
    result = recipe.create_result(master_fiberflat=rss, reduced_image=rss, reduced_rss=rss)
    result = recipe.run_qc(None, result)
    assert result.qc == QC.PARTIAL


def zero_order_image(amplitude, nfib=50, ncols=800):
    """2D flat with fibers every 7 rows, slightly tilted, and a diffuse band at column ~400"""
    nrows = 20 + 7 * nfib
    yy, xx = numpy.mgrid[0:nrows, 0:ncols].astype(float)
    traces = [(k + 1, numpy.polynomial.Polynomial([15.0 + 7 * k, 0.002])) for k in range(nfib)]
    spectrum = 1000.0 * (1 + 0.3 * xx / ncols)
    image = numpy.zeros((nrows, ncols))
    for _, poly in traces:
        image += spectrum * numpy.exp(-0.5 * ((yy - poly(xx)) / 1.2) ** 2)
    centre = 400 + 0.01 * (yy - nrows / 2)
    band = amplitude * 1000.0 * numpy.exp(-0.5 * ((xx - centre) / (14 / 2.3548)) ** 2)
    return image + band, traces, nrows


def test_find_zero_order():
    """A diffuse band at fixed columns, also between the fibers, is found with its centre and width"""
    image, traces, nrows = zero_order_image(0.1)
    bands = find_zero_order(image, traces)
    assert len(bands) == 1
    band = bands[0]
    assert abs(band["center"](nrows / 2) - 400) < 1
    assert abs(band["center"](nrows / 2 + 100) - 401) < 1
    assert 12 < band["fwhm"] < 16
    assert 0.07 < band["amplitude"] < 0.13
    # without the band, or with a band below the threshold, nothing is found
    assert find_zero_order(zero_order_image(0.0)[0], traces) == []
    assert find_zero_order(zero_order_image(0.01)[0], traces) == []


def test_remove_zero_order():
    """Each fiber is interpolated across the band, at the column of the band in its row"""
    ncols = 800
    xx = numpy.arange(ncols)
    traces = [(k + 1, numpy.polynomial.Polynomial([15.0 + 7 * k, 0.002])) for k in range(50)]
    band = {"center": numpy.polynomial.Polynomial([400.0 - 1.85, 0.01]), "fwhm": 14.0, "rows": (15.0, 358.0)}
    spectrum = 1000.0 * (1 + 0.3 * xx / ncols)
    rss = numpy.zeros((623, ncols))
    for fibid, poly in traces:
        col = band["center"](poly(400.0))
        rss[fibid - 1] = spectrum * (1 + 0.3 * numpy.exp(-0.5 * ((xx - col) / (14 / 2.3548)) ** 2))
    remove_zero_order(rss, traces, [band])
    ratio = rss[:50] / spectrum
    assert numpy.allclose(ratio, 1.0, atol=3e-3)
    # the rows without traces are not changed
    assert numpy.all(rss[50:] == 0)


def test_trace_polynomials():
    """The valid traces, with the global offset of the TraceMap"""
    import numpy.polynomial.polynomial as nppol
    from megaradrp.products.tracemap import GeometricTrace, TraceMap

    tracemap = TraceMap()
    tracemap.contents = [
        GeometricTrace(fibid=1, boxid=1, start=6, stop=4090, fitparms=[100.0, 0.01]),
        GeometricTrace(fibid=2, boxid=1, start=6, stop=4090, fitparms=[]),
    ]
    tracemap.global_offset = numpy.polynomial.Polynomial([0.5])
    traces = _trace_polynomials(tracemap)
    assert [fibid for fibid, _ in traces] == [1]
    assert numpy.isclose(traces[0][1](1000.0), nppol.polyval(1000.0, [100.0, 0.01]) + 0.5)


def test_trace_polynomials_modelmap():
    """With a ModelMap, the centre of the profiles; fibers without model are left out"""
    from megaradrp.products.modelmap import GeometricModel, ModelMap

    modelmap = ModelMap()
    centre = numpy.polynomial.Polynomial([200.0, 0.02])
    modelmap.contents = [
        GeometricModel(fibid=3, boxid=1, start=6, stop=4090, model={"params": {"mean": centre}}),
        GeometricModel(fibid=4, boxid=1, start=6, stop=4090, model={}),
    ]
    traces = _trace_polynomials(modelmap)
    assert [fibid for fibid, _ in traces] == [3]
    assert numpy.isclose(traces[0][1](1000.0), 220.0)
