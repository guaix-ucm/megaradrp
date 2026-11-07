#
# Copyright 2011-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline: Wavelength  Calibration"""

from numina.array.wavecalib.arccalibration import SolutionArcCalibration
import numpy.polynomial.polynomial as nppol

from .structured import BaseStructuredCalibration

# Flags of a wavelength solution that make it invalid
INVALID_SOLUTION_FLAGS = frozenset(["high_residual", "deviant"])


class FiberSolutionArcCalibration:
    """Wavelength solution of a fiber

    `flags` are the results of the quality checks of the solution:

    - "high_residual": the residual of the fit is large compared with the pixel
    - "deviant": the wavelengths differ from those of the neighbouring fibers
    - "edge_deviation": they differ only near the borders of the detector
    - "not_refined": the solution was not refined with all the lines

    A solution with a flag in INVALID_SOLUTION_FLAGS is not valid.
    """

    def __init__(self, fibid, solution, flags=None):
        self.fibid = fibid
        self.solution = solution
        self.flags = list(flags) if flags is not None else []

    @property
    def valid(self):
        return INVALID_SOLUTION_FLAGS.isdisjoint(self.flags)

    def add_flag(self, flag):
        """Add a quality flag to the solution, once"""
        if flag not in self.flags:
            self.flags.append(flag)

    def __getstate__(self):
        return {"fibid": self.fibid, "solution": self.solution.__getstate__(), "flags": self.flags}

    def __setstate__(self, state):
        self.fibid = state["fibid"]
        new = SolutionArcCalibration.__new__(SolutionArcCalibration)
        new.__setstate__(state["solution"])
        self.solution = new
        # solutions saved before the quality flags
        self.flags = list(state.get("flags", []))


class WavelengthCalibration(BaseStructuredCalibration):
    """Wavelength Calibration Product"""

    __tags__ = ["insmode", "vph"]

    def __init__(self, instrument="unknown"):
        super(WavelengthCalibration, self).__init__(instrument)
        self.contents = []
        self.global_offset = nppol.Polynomial([0.0])

    def tag_names(self):
        return ["insmode", "vph"]

    def invalid_fibers(self):
        """Fibers with an invalid wavelength solution"""
        return [val.fibid for val in self.contents if not val.valid]

    def __getstate__(self):
        st = super(WavelengthCalibration, self).__getstate__()

        st["contents"] = [val.__getstate__() for val in self.contents]
        st["global_offset"] = self.global_offset.coef
        return st

    def __setstate__(self, state):
        super(WavelengthCalibration, self).__setstate__(state)

        self.global_offset = nppol.Polynomial(state.get("global_offset", [0.0]))

        self.contents = []
        # Handle dictionary
        if isinstance(state["contents"], dict):
            values = state["contents"].values()
        else:
            values = state["contents"]

        for val in values:
            new = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
            new.__setstate__(val)
            self.contents.append(new)
