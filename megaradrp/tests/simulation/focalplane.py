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


class FocalPlane(object):

    NAMES = {'UNSET': 3, 'LEFT': 2, 'RIGHT': 1, 'SET': 0}
    CODES = {3: 'UNSET', 2: 'LEFT', 1: 'RIGHT', 0: 'SET'}
    def __init__(self):

        self._cover_status = 3

        self._filter_s = lambda x: True
        self._cover_s = lambda x: 1.0

        self.fibid = []
        self.pos = []

    def set_cover(self, mode):
        """Cover in the focal plane."""

        if mode.upper() not in self.NAMES:
            raise ValueError('"%s" mode is not recognized')

        self._cover_status = self.NAMES[mode.upper()]

        self._cover_s = lambda x: 1.0

        if self._cover_status == 3:
            self._filter_s = lambda x: True
            self._cover_s = lambda x: 1.0
        elif self._cover_status == 2:
            self._filter_s = lambda pos: pos[0] >= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] > 0.0 else 0.5
        elif self._cover_status == 1:
            self._filter_s = lambda pos: pos[0] <= 0.0
            self._cover_s = lambda pos: 1.0 if pos[0] < 0.0 else 0.5
        else:
            self._filter_s = lambda pos: False
            self._cover_s = lambda pos: 0.0

    def connect_fibers(self, fibid, pos):
        self.fibid = fibid
        self.pos = pos

    def get_visible_fibers(self):
        p1 = [(fid, self._cover_s(pos)) for fid, pos in zip(self.fibid, self.pos) if self._filter_s(pos)]
        return p1

    def meta(self):
        return {'cover': self.CODES[self._cover_status]}
