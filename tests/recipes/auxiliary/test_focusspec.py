import types

from astropy.io import fits
import numpy
import numpy.polynomial.polynomial as nppol
import pytest
from numina.array.wavecalib.solutionarc import CrLinear, SolutionArcCalibration
from numina.core import DataFrame
from numina.exceptions import RecipeError

import megaradrp.recipes.auxiliary.focusspec as focusspec
from megaradrp.products.tracemap import GeometricTrace, TraceMap
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration, WavelengthCalibration
from megaradrp.recipes.auxiliary.focusspec import FocusSpectrographRecipe

FOCI = [1.0, 2.0, 3.0, 4.0, 5.0]
BEST_FOCUS = 3.2
# columns of the lines in the RSS, 0-based
LINES = numpy.arange(300, 3900, 400)


def fwhm(focus):
    return 3.0 + 0.5 * (focus - BEST_FOCUS) ** 2


def make_rss(focus):
    """A RSS with lines whose width depends on the focus"""
    rng = numpy.random.default_rng(int(focus * 10))
    x = numpy.arange(4096)
    sigma = fwhm(focus) / (2 * numpy.sqrt(2 * numpy.log(2)))
    row = 10.0 + sum(5000.0 * numpy.exp(-0.5 * ((x - line) / sigma) ** 2) for line in LINES)
    data = row + rng.normal(0, 1, size=(623, 4096))
    hdr = fits.Header({"FOCUS": focus})
    return fits.HDUList([fits.PrimaryHDU(data.astype("float32"), header=hdr)])


def make_frame(focus):
    return DataFrame(frame=fits.HDUList([fits.PrimaryHDU(numpy.zeros((2, 2)), header=fits.Header({"FOCUS": focus}))]))


def make_tracemap():
    tracemap = TraceMap(instrument="MEGARA")
    tracemap.total_fibers = 623
    for fibid in range(1, 624):
        tracemap.contents.append(GeometricTrace(fibid, 1, 1, 4096, fitparms=[10.0 + 6.5 * fibid, 0.0]))
    return tracemap


def make_wlcalib(fibids):
    wlcalib = WavelengthCalibration(instrument="MEGARA")
    for fibid in fibids:
        sol = SolutionArcCalibration([], [6000.0, 0.3], 0.0, CrLinear(0, 0, 0, 0, 0))
        wlcalib.contents.append(FiberSolutionArcCalibration(fibid, sol))
    return wlcalib


class Extractor:
    """The extraction returns the RSS of the focus of the image"""

    def __init__(self, tracemap, datamodel, offset):
        pass

    def __call__(self, img):
        return make_rss(img[0].header["FOCUS"])


@pytest.fixture
def run_recipe(monkeypatch):
    """Run the recipe with the reduction and the extraction replaced"""

    def run(foci=FOCI, failed=(), vph="LR-R", **kwds):
        def processing(frames, flow, method, errors):
            with frames[0].open() as img:
                focus = img[0].header["FOCUS"]
            if focus in failed:
                raise ValueError("cannot be processed")
            return fits.HDUList([fits.PrimaryHDU(numpy.zeros((2, 2)), header=fits.Header({"FOCUS": focus}))])

        recipe = FocusSpectrographRecipe()
        monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: None)
        monkeypatch.setattr(focusspec, "basic_processing_with_combination_frames", processing)
        monkeypatch.setattr(focusspec, "ApertureExtractor", Extractor)
        # two images of each focus
        frames = [make_frame(focus) for focus in foci for _ in range(2)]
        obresult = types.SimpleNamespace(frames=frames, configuration=None, tags={"vph": vph, "insmode": "LCB"})
        values = dict(
            obresult=obresult,
            master_apertures=make_tracemap(),
            master_wlcalib=make_wlcalib(range(1, 624)),
            extraction_offset=[0.0],
            nfibers=100,
            tsigma=50,
        )
        values.update(kwds)
        return recipe.run.__wrapped__(recipe, types.SimpleNamespace(**values))

    return run


@pytest.mark.slow
def test_focusspec(run_recipe):
    result = run_recipe()
    # one measure for each line in fibers 1, 101, ..., 601
    final = result.focus_table
    assert final.shape == (7 * len(LINES), 3)
    assert numpy.allclose(final[:, 2], BEST_FOCUS, atol=0.05)
    # x of the line in the RSS, y of the trace
    assert numpy.allclose(sorted(set(numpy.round(final[:, 0]))), LINES)
    assert set(final[:, 1]) == {10.0 + 6.5 * fibid for fibid in range(1, 624, 100)}
    image = result.focus_image.open()[0].data
    assert image.shape == (4112, 4096)
    assert numpy.allclose(image, BEST_FOCUS, atol=0.05)
    # the measures of the lines in each image, with their wavelengths
    focus_wl = result.focus_wavelength
    assert sorted(focus_wl) == FOCI
    x, y, width, wl = focus_wl[1.0][101][0]
    assert width == pytest.approx(fwhm(1.0), abs=0.1)
    assert wl == pytest.approx(6000.0 + 0.3 * (4096 - x))


@pytest.fixture
def no_image(monkeypatch):
    """The image of the focus is slow to compute, replace it"""
    monkeypatch.setattr(FocusSpectrographRecipe, "generate_image", lambda self, final: numpy.zeros((2, 2)))


def test_focusspec_vph_without_limit(run_recipe, no_image, caplog):
    caplog.set_level("INFO")
    result = run_recipe(vph="XX-X", nfibers=300)
    assert "flux limit not defined for XX-X, using 40000.00" in caplog.text
    assert numpy.allclose(result.focus_table[:, 2], BEST_FOCUS, atol=0.05)


def test_focusspec_one_focus(run_recipe):
    with pytest.raises(RecipeError, match="only 1 different focus"):
        run_recipe(foci=[1.0])


@pytest.mark.xfail(reason="the foci that cannot be processed are used in the fit", strict=True)
def test_focusspec_focus_not_processed(run_recipe, no_image):
    """The foci that cannot be processed are ignored"""
    result = run_recipe(failed=[5.0], nfibers=300)
    assert numpy.allclose(result.focus_table[:, 2], BEST_FOCUS, atol=0.05)


def test_focusspec_wlcalib_fibers(run_recipe, no_image, caplog):
    """Only the fibers with wavelength calibration are measured"""
    result = run_recipe(master_wlcalib=make_wlcalib([50, 150]), nfibers=1)
    assert set(result.focus_table[:, 1]) == {10.0 + 6.5 * 50, 10.0 + 6.5 * 150}


def test_reorder_and_fit():
    recipe = FocusSpectrographRecipe()
    line_fibers = {
        1: {0: {"basedata": {"coordinates": (10.0, 20.0)}, "centers": [fwhm(f) for f in FOCI]}},
        2: {
            0: {"basedata": {"coordinates": (30.0, 40.0)}, "centers": [(f - 2.0) ** 2 + 1 for f in FOCI]},
            1: {"basedata": {"coordinates": (50.0, 60.0)}, "centers": [(f - 4.0) ** 2 + 1 for f in FOCI]},
        },
    }
    final = recipe.reorder_and_fit(line_fibers, FOCI)
    expected = [[10.0, 20.0, BEST_FOCUS], [30.0, 40.0, 2.0], [50.0, 60.0, 4.0]]
    assert numpy.allclose(final, expected)


def test_filter_lines():
    """The lines are matched with those of the image of the central focus"""
    recipe = FocusSpectrographRecipe()
    measures = {
        1.0: {1: [(100.0, 10.0, 5.0), (500.0, 10.0, 6.0)], 2: [(300.0, 20.0, 5.0)]},
        2.0: {1: [(100.5, 10.0, 4.0), (500.0, 10.0, 4.0), (900.0, 10.0, 4.0)], 2: []},
        # a line too far from the reference
        3.0: {1: [(101.0, 10.0, 5.5), (503.0, 10.0, 5.0)], 2: [(300.0, 20.0, 5.0)]},
    }
    line_fibers = recipe.filter_lines(measures)
    # fiber 2 has no lines in the reference
    assert list(line_fibers) == [1]
    # only the line at 100 is in the three images
    assert list(line_fibers[1]) == [0]
    assert line_fibers[1][0]["basedata"]["coordinates"] == (100.5, 10.0)
    assert line_fibers[1][0]["centers"] == [5.0, 4.0, 5.5]


def test_generate_focus_wl(caplog):
    recipe = FocusSpectrographRecipe()
    measures = {1.0: {1: [(96.0, 10.0, 5.0)], 7: [(96.0, 10.0, 5.0)]}}
    result = recipe.generate_focus_wl(measures, make_wlcalib([1]))
    assert result == {1.0: {1: [[96.0, 10.0, 5.0, nppol.polyval(4000.0, [6000.0, 0.3])]], 7: []}}
    assert "Fiber 7 hasn't WL calibration" in caplog.text
