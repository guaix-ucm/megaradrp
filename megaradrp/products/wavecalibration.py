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


from numina.array.wavecalib.arccalibration import SolutionArcCalibration

from .structured import BaseStructuredCalibration


class WavelengthCalibration(BaseStructuredCalibration):
    def __init__(self, instrument='unknown'):
        super(WavelengthCalibration, self).__init__(instrument)
        self.contents = {}

    def __getstate__(self):
        st = super(WavelengthCalibration, self).__getstate__()
        st['contents'] = {key: val.__getstate__()
                        for (key, val) in self.contents.items()}
        return st

    def __setstate__(self, state):
        super(WavelengthCalibration, self).__setstate__(state)
        self.contents = {}
        for (key, val) in state['contents'].items():
            new = SolutionArcCalibration.__new__(SolutionArcCalibration)
            new.__setstate__(val)
            self.contents[key] = new
        return self
