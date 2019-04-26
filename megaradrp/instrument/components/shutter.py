#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import six

from numina.instrument.components.wheel import Wheel
from numina.instrument.simulation.optics import Open, Stop, Filter


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
