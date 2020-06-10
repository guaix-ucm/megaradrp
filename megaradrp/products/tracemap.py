#
# Copyright 2011-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""

import numpy.polynomial.polynomial as nppol

from megaradrp.datatype import MegaraDataType
from .structured import BaseStructuredCalibration
from .aperture import GeometricAperture
from .traces import to_ds9_reg as to_ds9_reg_function


class GeometricTrace(GeometricAperture):
    """Representation of fiber trace on the image"""
    def __init__(self, fibid, boxid, start, stop, fitparms=None):
        super(GeometricTrace, self).__init__(fibid, boxid, start, stop)
        self.fitparms = fitparms if fitparms is not None else []
        self.polynomial = None
        # Update polynomial
        self._set_polynomial(fitparms)

    @property
    def valid(self):
        return self.is_valid()

    def is_valid(self):
        if self.fitparms:
            return True
        else:
            return False

    def __getstate__(self):
        state = super(GeometricTrace, self).__getstate__()
        del state['polynomial']
        return state

    def __setstate__(self, state):
        super(GeometricTrace, self).__setstate__(state)

        self._set_polynomial(state['fitparms'])

    def _set_polynomial(self, fitparms):
        if fitparms:
            self.polynomial = nppol.Polynomial(self.fitparms)
        else:
            self.polynomial = nppol.Polynomial([0.0])

    def aper_center(self):
        return self.polynomial


class TraceMap(BaseStructuredCalibration):
    """Trace map calibration product"""
    DATATYPE = MegaraDataType.TRACE_MAP
    __tags__ = ['insmode', 'vph']

    def __init__(self, instrument='MEGARA'):
        super(TraceMap, self).__init__(instrument)
        self.contents = []
        self.boxes_positions = []
        self.global_offset = nppol.Polynomial([0.0])
        self.ref_column = 2000
        self.expected_range = [4, 4092]
        #

    def __getstate__(self):
        st = super(TraceMap, self).__getstate__()
        st['contents'] = [t.__getstate__() for t in self.contents]
        st['boxes_positions'] = self.boxes_positions
        st['global_offset'] = self.global_offset.coef
        st['ref_column'] = self.ref_column
        st['expected_range'] = self.expected_range
        return st

    def __setstate__(self, state):
        super(TraceMap, self).__setstate__(state)
        self.contents = [GeometricTrace(**trace) for trace in state['contents']]
        # fibers in missing fibers and error_fitting are invalid
        self.boxes_positions = state.get('boxes_positions', [])
        self.global_offset = nppol.Polynomial(state.get('global_offset', [0.0]))
        self.ref_column = state.get('ref_column', 2000)
        self.expected_range = state.get('expected_range', [4, 4092])
        return self

    def to_ds9_reg(self, ds9reg, rawimage=False, numpix=100, fibid_at=0):
        """Transform fiber traces to ds9-region format.

        Parameters
        ----------
        ds9reg : BinaryIO
            Handle to output file name in ds9-region format.
        rawimage : bool
            If True the traces must be generated to be overplotted on
            raw FITS images.
        numpix : int
            Number of abscissae per fiber trace.
        fibid_at : int
            Abscissae where the fibid is shown (default=0 -> not shown).

        """
        return to_ds9_reg_function(self, ds9reg,
                                   rawimage=rawimage,
                                   numpix=numpix,
                                   fibid_at=fibid_at
                                   )
