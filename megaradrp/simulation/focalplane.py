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

import numpy


class FocalPlane(object):

    def __init__(self, cover):

        self.fiberset = {}
        self.cover = cover

    def set_cover(self, mode):
        """Cover in the focal plane."""

        self.cover.set(mode)

    def connect_lcb(self, lcb):
        self.fiberset['LCB'] = lcb

    def connect_fibermos(self, mos):
        self.fiberset['MOS'] = mos

    def get_visible_fibers(self, bundle):
        fibid, allpos = self.fiberset[bundle.name].fibers_in_focal_plane()

        p1 = self.cover.visible_fibers(fibid, allpos)

        return numpy.array(p1, dtype=[('fibid', 'i4'),('x', 'f4'), ('y', 'f4'), ('cover', 'f4')])

    def get_all_fibers(self, bundle):

        fibid, all_pos = self.fiberset[bundle.name].fibers_in_focal_plane()

        p2 = numpy.empty((len(fibid,)), dtype=[('fibid', 'i4'),('x', 'f4'), ('y', 'f4'), ('cover', 'f4')])
        p2['fibid'] = fibid
        p2['x'] = all_pos[:,0]
        p2['y'] = all_pos[:,1]
        p2['cover'] = 1.0

        return p2

    def config_info(self):
        return {'cover': self.cover.config_info()}
