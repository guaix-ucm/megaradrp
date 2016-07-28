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

import six
import numpy

from .wheel import Wheel


#
# These are optical elements, they should be implemented
# elsewhere


class Stop(object):
    def __init__(self, name):
        self.name = name

    def transmission(self, wl):
        return numpy.zeros_like(wl)


class Open(object):
    def __init__(self, name):
        self.name = name

    def transmission(self, wl):
        return numpy.ones_like(wl)


class Filter(object):
    def __init__(self, name, transmission=None):
        self.name = name

    def transmission(self, wl):
        # FIXME: implement this with a proper
        # transmission
        return numpy.ones_like(wl)


class MegaraShutter(Wheel):
    def __init__(self, parent=None):
        super(MegaraShutter, self).__init__(capacity=3, name='Shutter', parent=parent)
        self.put_in_pos(Stop(name='STOP'), 0) # FIXME
        self.put_in_pos(Open(name='OPEN'), 1) # FIXME
        # sorting order filter
        self.put_in_pos(Filter(transmission=None, name='FILTER'), 2) # FIXME
        self.move_to(1) # Open by default

        # MEGARA shutter has three positions:
        # open
        # closed
        # filter (sort order)

    def configure(self, value):
        # Let's see what is value:
        # a string
        if isinstance(value, six.string_types):
            val = value.lower()
            if val == 'open':
                val = 1
            elif val == 'closed':
                val = 0
            elif val == 'filter':
                val = 2
            else:
                raise ValueError('Not allowed value %s', value)
        elif isinstance(value, (int, long)):
            val = value
        else:
            raise TypeError('Not allowed type %s', type(value))

        # Move to value
        self.move_to(val)

    def open(self):
        self.move_to(1)

    def filter(self):
        self.move_to(2)

    def close(self):
        self.move_to(0)
