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
        self.state = state['fibid']
        new = SolutionArcCalibration.__new__(SolutionArcCalibration)
        new.__setstate__(state['solution'])
        self.solution = new


class WavelengthCalibration(BaseStructuredCalibration):
    def __init__(self, instrument='unknown'):
        super(WavelengthCalibration, self).__init__(instrument)
        self.total_fibers = 623
        self.missing_fibers = []
        self.error_fitting = []
        self.contents = {}

    def __getstate__(self):
        st = super(WavelengthCalibration, self).__getstate__()
        st['total_fibers'] = self.total_fibers
        st['missing_fibers'] = self.missing_fibers
        st['error_fitting'] = self.error_fitting

        st['contents'] = {key: val.__getstate__()
                        for (key, val) in self.contents.items()}
        return st

    def __setstate__(self, state):
        super(WavelengthCalibration, self).__setstate__(state)

        for key in state:
            if key not in ['contents']:
                setattr(self, key, state[key])

        self.contents = {}
        for (key, val) in state['contents'].items():
            new = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
            new.__setstate__(val)
            self.contents[key] = new

    @property
    def calibid(self):
        return 'uuid:{}'.format(self.uuid)

    @classmethod
    def _datatype_dump(cls, obj, where):
        filename = where.destination + '.json'

        with open(filename, 'w') as fd:
            json.dump(obj.__getstate__(), fd, indent=2)

        return filename

    @classmethod
    def _datatype_load(cls, obj):
        try:
            with open(obj, 'r') as fd:
                state = json.load(fd)
        except IOError as e:
            raise e

        result = cls.__new__(cls)
        result.__setstate__(state=state)
        return result
