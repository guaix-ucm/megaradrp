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

"""Products of the Megara Pipeline: Wavelength  Calibration"""

import json

from numina.array.wavecalib.arccalibration import SolutionArcCalibration
import numpy.polynomial.polynomial as nppol

from .structured import BaseStructuredCalibration


class FiberSolutionArcCalibration(object):
    def __init__(self, fibid, solution):
        self.fibid = fibid
        self.solution = solution

    def __getstate__(self):
        return {'fibid': self.fibid,
                'solution': self.solution.__getstate__()
        }

    def __setstate__(self, state):
        self.fibid = state['fibid']
        new = SolutionArcCalibration.__new__(SolutionArcCalibration)
        new.__setstate__(state['solution'])
        self.solution = new


class WavelengthCalibration(BaseStructuredCalibration):
    """Wavelength Calibration Product
    """
    def __init__(self, instrument='unknown'):
        super(WavelengthCalibration, self).__init__(instrument)
        self.contents = []
        self.global_offset = nppol.Polynomial([0.0])

    def __getstate__(self):
        st = super(WavelengthCalibration, self).__getstate__()

        st['contents'] = [val.__getstate__() for val in self.contents]
        st['global_offset'] = self.global_offset.coef
        return st

    def __setstate__(self, state):
        super(WavelengthCalibration, self).__setstate__(state)

        self.global_offset = nppol.Polynomial(state.get('global_offset', [0.0]))

        self.contents = []
        # Handle dictionary
        if isinstance(state['contents'], dict):
            values = state['contents'].values()
        else:
            values = state['contents']

        for val in values:
            new = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
            new.__setstate__(val)
            self.contents.append(new)
