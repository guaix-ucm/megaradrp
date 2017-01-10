#
# Copyright 2016 Universidad Complutense de Madrid
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
import uuid

import numina.core.types
import numina.core.products


class BaseStructuredCalibration(numina.core.products.DataProductTag,
                                numina.core.types.AutoDataType):
    """Base class for structured calibration data

    Parameters
    ----------

    instrument: str
        Instrument name

    Attributes
    ----------
    tags: dict
        dictionary of selection fields
    uuid: str
       UUID of the result

    """
    def __init__(self, instrument='unknown'):
        super(BaseStructuredCalibration, self).__init__()
        self.instrument = instrument
        self.tags = {}
        self.uuid = uuid.uuid1().hex
        self.total_fibers = 0
        self.missing_fibers = []
        self.error_fitting = []

    @property
    def calibid(self):
        return 'uuid:{}'.format(self.uuid)

    @property
    def default(self):
        return None

    def __getstate__(self):
        st = {}
        for key in ['instrument', 'tags', 'uuid']:
            st[key] = self.__dict__[key]
        st['total_fibers'] = self.total_fibers
        st['missing_fibers'] = self.missing_fibers
        st['error_fitting'] = self.error_fitting
        return st

    def __setstate__(self, state):
        self.instrument = state['instrument']
        self.tags = state['tags']
        self.uuid = state['uuid']

        for key in state:
            if key not in ['contents']:
                setattr(self, key, state[key])

    def __str__(self):
        sclass = type(self).__name__
        return "{}(instrument={}, uuid={})".format(sclass, self.instrument, self.uuid)

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
