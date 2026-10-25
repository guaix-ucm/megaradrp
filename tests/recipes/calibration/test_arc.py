#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Tests for the arc mode recipe module."""

import logging
import pathlib
import types

from astropy.io import fits
import numpy
import numpy.polynomial.polynomial as nppol
import pytest
from numina.array.wavecalib.solutionarc import CrLinear, SolutionArcCalibration
from numina.core import DataFrame
from numina.user.cli import main
from numina.util.node import Node

import megaradrp.recipes.calibration.arc as arcmod
from megaradrp.loader import load_drp
from megaradrp.products.tracemap import GeometricTrace, TraceMap
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration, WavelengthCalibration
from megaradrp.recipes.calibration.arc import ArcCalibrationRecipe
from megaradrp.testing.common import make_obresult, raw_header

BASE_URL = "https://guaix.fis.ucm.es/~spr/megara_test/"


def run_recipe():

    main(["run", "obsrun.yaml", "-r", "control.yaml"])


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_arc_calibration_set0(drpmocker):

    drpmocker.add_drp("MEGARA", load_drp)

    run_recipe()


# The wavelength calibration with a synthetic arc


@pytest.fixture(autouse=True)
def no_debug_plots(caplog):
    """find_fxpeaks plots each spectrum with the level DEBUG, that numina run leaves set"""
    caplog.set_level(logging.INFO, logger="numina.array.wavecalib")


DATA = pathlib.Path(__file__).parent / "data"


def wl_true(channel, fibid):
    """The wavelength of the channel (1-based pixel) of the fiber"""
    return 6030.0 + 0.002 * fibid + 0.31 * (channel - 1) - 2e-6 * (channel - 1) ** 2


# fibers with a trace, the others are missing
VALID = list(range(10, 130, 10))


def load_catalog(name):
    """Wavelength, flag and intensity of the lines"""
    return numpy.loadtxt(DATA / name, usecols=(0, 1, 2))


def arc_rss(catalog):
    """RSS of an arc with the lines of the catalog, with their intensities"""
    x = numpy.arange(4096)
    rss = numpy.zeros((623, 4096))
    for fibid in VALID:
        wl = wl_true(x + 1, fibid)
        for wavelength, _, intensity in catalog:
            center = numpy.interp(wavelength, wl, x)
            rss[fibid - 1] += intensity * numpy.exp(-0.5 * ((x - center) / 1.5) ** 2)
        rss[fibid - 1] += 5.0
    return rss


@pytest.fixture(scope="module")
def thne():
    """All the lines have flag 1"""
    catalog = load_catalog("LR-R_ThNe.lis")
    return catalog, arc_rss(catalog)


@pytest.fixture(scope="module")
def thar():
    """The lines with flag 0 are only used to refine the calibration"""
    catalog = load_catalog("LR-R_ThAr.lis")
    return catalog, arc_rss(catalog)


def make_tracemap(valid=VALID):
    tracemap = TraceMap(instrument="MEGARA")
    tracemap.total_fibers = 623
    for fibid in range(1, 624):
        fitparms = [10.0 + 6.5 * fibid, 0.0] if fibid in valid else []
        tracemap.contents.append(GeometricTrace(fibid, 1, 1, 4096, fitparms=fitparms))
    return tracemap


@pytest.fixture
def recipe(monkeypatch):
    """The recipe, the image of the FWHM is slow to compute"""
    recipe = ArcCalibrationRecipe()
    monkeypatch.setattr(recipe, "generate_fwhm_image", lambda solutions: numpy.zeros((2, 2)))
    return recipe


def wl_error(solution, fibid):
    channels = numpy.array([1.0, 2000.0, 4096.0])
    return numpy.abs(nppol.polyval(channels, solution.coeff) - wl_true(channels, fibid)).max()


def test_calibrate_wl(thne, recipe):
    """Without lines with flag 0, there is no refinement"""
    catalog, rss = thne
    # fiber 5 has a trace but no lines
    initial, refined, fwhm_image = recipe.calibrate_wl(rss, catalog, [3], make_tracemap(VALID + [5]), [0])
    assert refined is None
    assert initial.error_fitting == [5]
    assert len(initial.missing_fibers) == 623 - len(VALID) - 1
    assert [sol.fibid for sol in initial.contents] == VALID
    for fibsol in initial.contents:
        assert wl_error(fibsol.solution, fibsol.fibid) < 1e-3
        # the position of the line in the trace, and its FWHM
        for feature in fibsol.solution.features:
            assert feature.ypos == 10.0 + 6.5 * fibsol.fibid
            assert feature.fwhm > 0
    assert fwhm_image[0].data.shape == (2, 2)


@pytest.mark.xfail(reason="the FWHM is measured around the 1-based channel of the line", strict=True)
def test_calibrate_wl_fwhm(thne, recipe):
    """The FWHM of the lines is measured"""
    catalog, rss = thne
    initial, _, _ = recipe.calibrate_wl(rss, catalog, [3], make_tracemap([10]), [0])
    features = [feature for feature in initial.contents[0].solution.features if 100 < feature.xpos < 4000]
    # the FWHM is measured around the pixel of the peak, it is overestimated up to 0.2
    # pixels when the line is in the border of a pixel
    assert numpy.median([feature.fwhm for feature in features]) == pytest.approx(1.5 * 2.3548, abs=0.15)


def test_calibrate_wl_old_format(thne, recipe):
    """The catalog of two columns has no flags, there is no refinement"""
    catalog, rss = thne
    initial, refined, _ = recipe.calibrate_wl(rss, catalog[:, [0, 2]], [3], make_tracemap(), [0])
    assert refined is None
    assert all(wl_error(sol.solution, sol.fibid) < 1e-3 for sol in initial.contents)


def test_calibrate_wl_refined(thar, recipe):
    """The lines with flag 0 are used to refine the calibration"""
    catalog, rss = thar
    initial, refined, _ = recipe.calibrate_wl(rss, catalog, [3, 3], make_tracemap(), [30])
    assert [sol.fibid for sol in initial.contents] == VALID
    assert [sol.fibid for sol in refined.contents] == VALID
    assert len(refined.missing_fibers) == 623 - len(VALID)
    assert refined.total_fibers == 623
    nflag1 = int(catalog[:, 1].sum())
    for fibsol in refined.contents:
        assert wl_error(fibsol.solution, fibsol.fibid) < 1e-3
        # more lines than those with flag 1
        assert fibsol.solution.npoints_eff > nflag1


@pytest.mark.parametrize(
    "ncols, poldeg, message",
    [
        (2, [3, 4, 5], "Unexpected list length"),
        (2, [4, 3], "Unexpected poldeg_initial"),
        (4, [3], "expected number of columns"),
    ],
)
def test_calibrate_wl_errors(thne, recipe, ncols, poldeg, message):
    catalog, rss = thne
    catalog = numpy.column_stack([catalog[:, 0]] + [catalog[:, 1]] * (ncols - 1))
    with pytest.raises(ValueError, match=message):
        recipe.calibrate_wl(rss, catalog, poldeg, make_tracemap(), [0])


def test_model_coeff_vs_fiber_degree():
    wlcalib = WavelengthCalibration(instrument="MEGARA")
    sol = SolutionArcCalibration([], [6000.0, 0.3], 0.0, CrLinear(0, 0, 0, 0, 0))
    wlcalib.contents.append(FiberSolutionArcCalibration(1, sol))
    with pytest.raises(ValueError, match="Unexpected number of polynomial coefficients"):
        ArcCalibrationRecipe().model_coeff_vs_fiber(wlcalib, 3)


@pytest.mark.slow
def test_generate_fwhm_image():
    """The FWHM of the nearest line, of the fibers multiple of 10"""
    solutions = []
    for fibid, (xpos, ypos, fwhm) in [
        (10, (100.0, 200.0, 3.0)),
        (11, (3000.0, 3000.0, 9.0)),
        (20, (3000.0, 4000.0, 4.0)),
    ]:
        sol = SolutionArcCalibration([], [6000.0, 0.3], 0.0, CrLinear(0, 0, 0, 0, 0))
        sol.features = [types.SimpleNamespace(xpos=xpos, ypos=ypos, fwhm=fwhm)]
        solutions.append(FiberSolutionArcCalibration(fibid, sol))
    image = ArcCalibrationRecipe().generate_fwhm_image(solutions)
    assert image.shape == (4112, 4096)
    assert image[0, 0] == 3.0
    assert image[3000, 3000] == 4.0
    assert set(numpy.unique(image)) == {3.0, 4.0}


class Extractor(Node):
    """The extraction returns the RSS of the arc, before the flip"""

    rss = None

    def __init__(self, tracemap, datamodel, offset):
        super().__init__()

    def run(self, img):
        return fits.HDUList([fits.PrimaryHDU(self.rss[:, ::-1].astype("float32"), header=img[0].header)])


@pytest.fixture
def run_arc(monkeypatch, thar):
    """Run the recipe, with the reduction and the extraction replaced"""

    def run(recipe, vph="LR-R", **kwds):
        catalog, rss = thar
        reduced = fits.HDUList([fits.PrimaryHDU(numpy.zeros((10, 10)))])
        monkeypatch.setattr(recipe, "init_filters", lambda rinput, conf: None)
        monkeypatch.setattr(
            arcmod, "basic_processing_with_combination", lambda rinput, flow, method, method_kwargs: reduced
        )
        monkeypatch.setattr(Extractor, "rss", rss)
        monkeypatch.setattr(arcmod, "ApertureExtractor", Extractor)
        header = raw_header(vph=vph, insmode="LCB")
        obresult = make_obresult("MegaraArcCalibration", header=header)
        obresult.tags = {"vph": vph, "insmode": "LCB"}
        raw = fits.PrimaryHDU(header=header)
        raw.header["UUID"] = "11111111-2222-3333-4444-555555555555"
        obresult.frames = [DataFrame(frame=fits.HDUList([raw]))]
        values = dict(
            # the bias is not used, the reduction is replaced
            master_bias=DataFrame(frame=fits.HDUList([fits.PrimaryHDU()])),
            master_apertures=make_tracemap(),
            lines_catalog=catalog,
            polynomial_degree=[3, 3],
            nlines=[30],
        )
        values.update(kwds)
        rinput = recipe.create_input(obresult=obresult, **values)
        return recipe.run(rinput)

    return run


def test_arc_run(run_arc, recipe):
    """The RSS is flipped, the calibration is the refined one with the metadata of the initial"""
    result = run_arc(recipe)
    wlcalib = result.master_wlcalib
    assert [sol.fibid for sol in wlcalib.contents] == VALID
    assert all(wl_error(sol.solution, sol.fibid) < 1e-3 for sol in wlcalib.contents)
    assert wlcalib.tags == {"vph": "LR-R", "insmode": "LCB"}
    assert wlcalib.meta_info["origin"]["insconf_uuid"] == raw_header()["INSCONF"]
    rss = result.reduced_rss.open()[0].data
    assert numpy.allclose(rss[VALID[0] - 1], Extractor.rss[VALID[0] - 1])


def test_arc_run_vph_without_threshold(run_arc, recipe, caplog):
    """Without refinement, the result is the initial calibration"""
    caplog.set_level("INFO")
    result = run_arc(recipe, vph="XX-X", polynomial_degree=[3])
    assert "rel threshold not defined for XX-X" in caplog.text
    assert len(result.master_wlcalib.contents) == len(VALID)


@pytest.mark.slow
# the lines not identified have funcost=inf
@pytest.mark.filterwarnings("ignore:Converting funcost=inf:RuntimeWarning")
def test_arc_run_intermediate(run_arc, tmp_path, monkeypatch):
    """The plots of the calibration are saved"""
    monkeypatch.chdir(tmp_path)
    recipe = ArcCalibrationRecipe(intermediate_results=True)
    monkeypatch.setattr(recipe, "generate_fwhm_image", lambda solutions: numpy.zeros((2, 2)))
    run_arc(recipe, store_pdf_with_refined_fits=1)
    names = {path.name for path in tmp_path.iterdir()}
    expected = {"wavecal_iter1.pdf", "wavecal_iter2.pdf", "wavecal_refine_iter1.pdf", "initial_master_wlcalib.json"}
    assert expected <= names
    assert {path.name for path in (tmp_path / "refined_wavecal").iterdir()} == {f"{fibid:03d}.pdf" for fibid in VALID}
