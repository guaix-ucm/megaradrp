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

import uuid
import yaml

import numina.core.types
import numina.core.products


class GeometricTrace(object):
    def __init__(self, fibid, boxid, start, stop, fitparms=None):
        self.fibid = fibid
        self.boxid = boxid
        self.start = start
        self.stop = stop
        self.fitparms = fitparms if fitparms is not None else []

    def __getstate__(self):
        return self.__dict__


class TraceMap(numina.core.types.AutoDataType):
    def __init__(self, instrument='unknown'):
        super(TraceMap, self).__init__()

        self.instrument = instrument
        self.tags = {}
        self.uuid = uuid.uuid1().hex
        self.tracelist = []

    def __getstate__(self):
        st = {}
        st['tracelist'] = [t.__getstate__() for t in self.tracelist]
        for key in ['instrument', 'tags', 'uuid']:
            st[key] = self.__dict__[key]
        return st

    def __setstate__(self, state):
        self.instrument = state['instrument']
        self.tags = state['tags']
        self.uuid = state['uuid']
        self.tracelist = [GeometricTrace(**trace) for trace in state['tracelist']]

        return self

    @classmethod
    def _datatype_dump(cls, obj, where):
        filename = where.destination + '.yaml'

        with open(filename, 'w') as fd:
            yaml.dump(obj.__getstate__(), fd)

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
