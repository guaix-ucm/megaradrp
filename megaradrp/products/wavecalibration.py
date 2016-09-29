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


import uuid

import yaml
import json

import numina.core.types
import numina.core.products
from numina.array.wavecalib.arccalibration import SolutionArcCalibration

class BaseStructuredCalibration(numina.core.products.DataProductTag,
                                numina.core.types.AutoDataType):
    def __init__(self, instrument='unknown'):
        super(BaseStructuredCalibration, self).__init__()
        self.instrument = instrument
        self.tags = {}
        self.uuid = uuid.uuid1().hex

    def __getstate__(self):
        st = {}
        for key in ['instrument', 'tags', 'uuid']:
            st[key] = self.__dict__[key]
        return st


class WavelengthCalibration(BaseStructuredCalibration):
    def __init__(self, instrument='unknown'):
        super(WavelengthCalibration, self).__init__()
        self.contents = {}

    def __getstate__(self):
        st = {}
        st['contents'] = {key: val.__getstate__()
                        for (key, val) in self.contents.items()}
        st.update(BaseStructuredCalibration.__getstate__(self))
        return st

    def __setstate__(self, state):
        self.instrument = state['instrument']
        self.tags = state['tags']
        self.uuid = state['uuid']
        self.contents = {key: SolutionArcCalibration(**val)
                              for (key, val) in state['contents'].items()}
        return self

    @classmethod
    def _datatype_dump(cls, obj, where):
        filename = where.destination + '.yaml'

        with open('xxx.txt', 'w') as fd:
            fd.write(str(obj.__getstate__()))
        with open(filename, 'w') as fd:
            json.dump(obj.__getstate__(), fd, indent=4, sort_keys=True)

        return filename

    @classmethod
    def _datatype_load(cls, obj):

        try:
            with open(obj, 'r') as fd:
                state = yaml.load(fd)
        except IOError as e:
            raise e

        result = cls.__new__(cls)
        result.__setstate__(state=state)
        return result

    @property
    def default(self):
        return None
