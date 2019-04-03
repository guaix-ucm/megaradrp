#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy as np

from numina.instrument.hwdevice import HWDevice
from numina.instrument.components.signal import Signal


class HemiCover(HWDevice):
    STATES = [0, 1]
    NAMES = {'PARKED': 1, 'INPLACE': 0}
    CODES = {1: 'PARKED', 0: 'INPLACE'}

    def __init__(self, name, parent=None):
        self._pos = 0 # Closed
        super(HemiCover, self).__init__(name, parent=parent)

        self.changed = Signal()
        self.opened = Signal()
        self.closed = Signal()

    def set(self, pos):
        if pos not in self.STATES:
            raise ValueError('%d is not a valid state' % pos)
        if pos != self._pos: # We have to move
            self.flip()

    def open(self):
        if self._pos == 0:
            self._pos = 1
            self.changed.emit(self._pos)
            self.opened.emit()

    def close(self):
        if self._pos == 1:
            self._pos = 0
            self.changed.emit(self._pos)
            self.closed.emit()

    def flip(self):
        self._pos += 1
        self._pos %= 2
        self.changed.emit(self._pos)
        if self._pos == 1:
            self.opened.emit()
        else:
            self.closed.emit()

    @property
    def position(self):
        return self._pos

    @property
    def label(self):
        return self.CODES[self._pos]


class FullCover(HWDevice):
    # STATES = 00, 01, 10, 11 
    # FULL CLOSED, CLOSED LEFT, CLOSED RIGHT, FULL OPEN
    STATES = [0, 1, 2, 3]
    NAMES = {'UNSET': 3, 'LEFT': 2, 'RIGHT': 1, 'SET': 0}
    CODES = {3: 'UNSET', 2: 'LEFT', 1: 'RIGHT', 0: 'SET'}

    def __init__(self, name=None, parent=None):
        super(FullCover, self).__init__(name, parent=parent)
        self.left = HemiCover(name='Left', parent=self)
        self.right = HemiCover(name='Right', parent=self)

        self.changed_left = self.left.changed
        self.opened_left = self.left.opened
        self.closed_left = self.left.closed

        self.changed_right = self.right.changed
        self.opened_right = self.right.opened
        self.closed_right = self.right.closed

    def flip(self):
        self.left.flip()
        self.right.flip()
        
    def open(self):
        self.left.open()
        self.right.open()
        
    def close(self):
        self.left.close()
        self.right.close()

    def set(self, pos):
        if pos not in self.STATES:
            raise ValueError('%d is not a valid state' % pos)

        l_pos = pos // 2
        r_pos = pos % 2
        self.left.set(l_pos)
        self.right.set(r_pos)

    def set_mode(self, mode):
        pos = self.NAMES[mode]
        self.set(pos)

    @property
    def position(self):
        return 2 * self.left.position + self.right.position

    @position.setter
    def position(self, pos):
        self.set(pos)

    @property
    def label(self):
        return self.CODES[self.position]

    @label.setter
    def label(self, value):
        self.set_mode(value)


class MegaraCover(FullCover):
    """MEGARA Cover"""
    NAMES = {'UNSET': 3, 'LEFT': 2, 'RIGHT': 1, 'SET': 0}
    CODES = {3: 'UNSET', 2: 'LEFT', 1: 'RIGHT', 0: 'SET'}
    HCODES = {1: 'PARKED', 0: 'INPLACE'}
  
    VALS = {3: lambda pos: np.ones(pos[:,0].shape, dtype='bool'),
            2: lambda pos: pos[:,0]<-0.7,
            1: lambda pos: pos[:,0]>0.7,
            0: lambda pos: np.zeros(pos[:,0].shape, dtype='bool')}

    def __init__(self, parent=None):
        self.mode = 'UNSET'
        super(MegaraCover, self).__init__(name='Cover', parent=parent)
        self.set_mode(self.mode)

        self._filter_s = lambda x: True
        self._cover_s = lambda x: 1.0

    def set_mode(self, mode):
        """Cover in the focal plane."""
        if mode.upper() not in self.NAMES:
            raise ValueError('"%s" mode is not recognized' % mode)

        pos = self.NAMES[mode.upper()]
        super(MegaraCover, self).set(pos)

        assert pos == self.position

        self._cover_s = lambda x: 1.0

        if self.position == 3:
            self._filter_s = lambda x: True
            self._cover_s = lambda x: 1.0
        elif self.position == 2:
            self._filter_s = lambda pos: pos[0] >= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] > 0.0 else 0.5
        elif self.position == 1:
            self._filter_s = lambda pos: pos[0] <= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] < 0.0 else 0.5
        else:
            self._filter_s = lambda pos: False
            self._cover_s = lambda pos: 0.0

    def visible_fibers(self, fibid, allpos):
        p1 = [(fid, pos[0], pos[1], self._cover_s(pos)) for fid, pos in zip(fibid, allpos) if self._filter_s(pos)]
        return p1

    def __call__(self, fiber_positions):
        return self.VALS[self.position](fiber_positions)
