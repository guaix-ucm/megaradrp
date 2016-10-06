#
# Copyright 2011-2016 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

"""Products of the Megara Pipeline"""

import numpy.polynomial.polynomial as nppol

from .structured import BaseStructuredCalibration


class GeometricTrace(object):
    def __init__(self, fibid, boxid, start, stop, fitparms=None):
        self.fibid = fibid
        self.boxid = boxid
        self.start = start
        self.stop = stop
        self.fitparms = fitparms if fitparms is not None else []
        self.polynomial = None
        # Update polynomial
        self._set_polynomial(fitparms)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['polynomial']
        return state

    def __setstate__(self, state):
        self.__dict__  = state
        self._set_polynomial(state['fitparms'])

    def _set_polynomial(self, fitparms):
        if fitparms:
            self.polynomial = nppol.Polynomial(self.fitparms)
        else:
            self.polynomial = nppol.Polynomial([0.0])



class TraceMap(BaseStructuredCalibration):
    def __init__(self, instrument='unknown'):
        super(TraceMap, self).__init__(instrument)
        self.contents = []

    def __getstate__(self):
        st = super(TraceMap, self).__getstate__()
        st['contents'] = [t.__getstate__() for t in self.contents]
        return st

    def __setstate__(self, state):
        super(TraceMap, self).__setstate__(state)
        self.contents = [GeometricTrace(**trace) for trace in state['contents']]

        return self
