#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Quality checks of the wavelength solutions of the fibers

The wavelength solutions of neighbouring fibers are almost equal. A solution
with a large residual, or whose wavelengths differ from those of the
neighbouring fibers, has lines identified wrongly.
"""

import logging

import numpy
import numpy.polynomial.polynomial as nppol

import numina.types.qc as qc

_logger = logging.getLogger(__name__)


def flag_solutions(wlcalib, naxis1, refined=True, max_residual=0.5, max_deviation=4.0, nneigh=5, margin=0.1, ncols=9):
    """Add the quality flags to the wavelength solutions

    - "high_residual": the residual of the fit is above `max_residual` pixels
    - "deviant": in the interior of the detector (leaving `margin` times `naxis1`
      at each side), the wavelengths differ from the median of those of the
      `nneigh` fibers at each side by more than `max_deviation` pixels
    - "edge_deviation": they differ by more only at the borders of the detector,
      where the solutions extrapolate beyond the lines
    - "not_refined": the solutions are those of the first pass (`refined` is False)

    Parameters
    ----------
    wlcalib : WavelengthCalibration
        The solutions, changed in place
    naxis1 : int
        Number of pixels of the spectra
    refined : bool
        The solutions were refined with all the lines
    """
    items = sorted((item for item in wlcalib.contents if len(item.solution.coeff) > 1), key=lambda item: item.fibid)
    if not items:
        return
    xinner = numpy.linspace(1 + margin * naxis1, naxis1 - margin * naxis1, ncols)
    xedges = numpy.array([1.0, float(naxis1)])
    lam_inner = numpy.array([nppol.polyval(xinner, item.solution.coeff) for item in items])
    lam_edges = numpy.array([nppol.polyval(xedges, item.solution.coeff) for item in items])
    for k, item in enumerate(items):
        solution = item.solution
        cdelt = abs(solution.cr_linear.cdelt)
        if not refined:
            item.add_flag("not_refined")
        if cdelt == 0:
            continue
        if solution.residual_std / cdelt > max_residual:
            item.add_flag("high_residual")
        neighbours = [j for j in range(max(0, k - nneigh), min(len(items), k + nneigh + 1)) if j != k]
        if not neighbours:
            continue
        dev_inner = numpy.max(numpy.abs(lam_inner[k] - numpy.median(lam_inner[neighbours], axis=0))) / cdelt
        dev_edges = numpy.max(numpy.abs(lam_edges[k] - numpy.median(lam_edges[neighbours], axis=0))) / cdelt
        if dev_inner > max_deviation:
            item.add_flag("deviant")
        elif dev_edges > max_deviation:
            item.add_flag("edge_deviation")


def check_qc_wlcalib(wlcalib, nfibers, ignored_ids=(), logger=None):
    """Global quality of the wavelength calibration

    GOOD if every fiber with spectrum has a valid solution, PARTIAL if some do not,
    BAD if more than 5 % of the fibers do not. The fibers without spectrum
    (`missing_fibers`) do not count.
    """
    logger = logger or _logger
    failed = sorted(set(wlcalib.error_fitting) - set(ignored_ids))
    invalid = sorted(set(wlcalib.invalid_fibers()) - set(ignored_ids))
    counts = {}
    for item in wlcalib.contents:
        for flag in item.flags:
            counts.setdefault(flag, []).append(item.fibid)

    if failed:
        logger.warning("fibers without wavelength solution: %s", failed)
    if invalid:
        logger.warning("fibers with an invalid wavelength solution: %s", invalid)
    for flag, fibids in sorted(counts.items()):
        if flag == "not_refined":
            logger.info("the wavelength solutions are not refined")
        else:
            logger.info("solutions with flag %r: %d %s", flag, len(fibids), fibids)

    nbad = len(failed) + len(invalid)
    if nbad > 0.05 * nfibers:
        wlcalib.quality_control = qc.QC.BAD
    elif nbad > 0:
        wlcalib.quality_control = qc.QC.PARTIAL
    else:
        wlcalib.quality_control = qc.QC.GOOD
    return wlcalib.quality_control
