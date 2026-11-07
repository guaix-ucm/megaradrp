import logging

import numpy.polynomial.polynomial as nppol
import pytest
from numina.array.wavecalib.solutionarc import CrLinear, SolutionArcCalibration
import numina.types.qc as qc

from megaradrp.processing.arcqc import check_qc_wlcalib, flag_solutions
from megaradrp.products import WavelengthCalibration
from megaradrp.products.wavecalibration import FiberSolutionArcCalibration

NAXIS1 = 4096
CDELT = 0.37


def make_solution(coeff, residual=0.01):
    poly = nppol.Polynomial(coeff)
    crmin, crmax = poly(1), poly(NAXIS1)
    cr = CrLinear(crpix=1.0, crval=crmin, crmin=crmin, crmax=crmax, cdelt=(crmax - crmin) / (NAXIS1 - 1))
    return SolutionArcCalibration(features=[], coeff=list(coeff), residual_std=residual, cr_linear=cr)


def make_wlcalib(nfibers=60, special=None):
    special = special or {}
    wlcalib = WavelengthCalibration(instrument="MEGARA")
    for fibid in range(1, nfibers + 1):
        coeff, residual = special.get(fibid, ([7000.0 + 0.01 * fibid, CDELT], 0.01))
        wlcalib.contents.append(FiberSolutionArcCalibration(fibid, make_solution(coeff, residual)))
    return wlcalib


def edge_coeff(fibid, dev_edge_px=8.0):
    """A solution that departs from the others only near the borders of the detector"""
    a = dev_edge_px * CDELT / (NAXIS1 / 2) ** 4
    quartic = (nppol.Polynomial([-NAXIS1 / 2, 1.0]) ** 4).coef * a
    base = [7000.0 + 0.01 * fibid, CDELT, 0, 0, 0]
    return [b + q for b, q in zip(base, quartic)]


def test_flag_solutions():
    special = {
        30: ([7010.3, CDELT], 0.01),  # ~28 px away from its neighbours
        40: ([7000.4, CDELT], 0.5),  # residual of 1.35 px
        50: (edge_coeff(50), 0.01),
    }
    wlcalib = make_wlcalib(special=special)
    flag_solutions(wlcalib, NAXIS1)
    flags = {item.fibid: item.flags for item in wlcalib.contents if item.flags}
    assert flags == {30: ["deviant"], 40: ["high_residual"], 50: ["edge_deviation"]}
    assert sorted(wlcalib.invalid_fibers()) == [30, 40]


def test_flag_solutions_not_refined():
    wlcalib = make_wlcalib(nfibers=20)
    flag_solutions(wlcalib, NAXIS1, refined=False)
    assert all(item.flags == ["not_refined"] for item in wlcalib.contents)
    assert wlcalib.invalid_fibers() == []


@pytest.mark.parametrize(
    "special, failed, expected",
    [
        ({}, [], qc.QC.GOOD),
        ({30: ([7010.3, CDELT], 0.01)}, [], qc.QC.PARTIAL),
        ({}, [623], qc.QC.GOOD),
        ({}, list(range(100, 140)), qc.QC.BAD),
    ],
)
def test_check_qc_wlcalib(special, failed, expected, caplog):
    wlcalib = make_wlcalib(special=special)
    wlcalib.error_fitting = failed
    flag_solutions(wlcalib, NAXIS1)
    with caplog.at_level(logging.INFO):
        result = check_qc_wlcalib(wlcalib, nfibers=623, ignored_ids=[623])
    assert result == expected
    assert wlcalib.quality_control == expected


def test_solution_flags_state():
    """The flags are kept in the state; solutions saved before have no flags"""
    item = FiberSolutionArcCalibration(1, make_solution([7000.0, CDELT]), flags=["deviant"])
    state = item.__getstate__()
    new = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
    new.__setstate__(state)
    assert new.flags == ["deviant"]
    assert not new.valid
    del state["flags"]
    old = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
    old.__setstate__(state)
    assert old.flags == []
    assert old.valid
