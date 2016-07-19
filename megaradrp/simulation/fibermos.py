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

import inspect

from .device import HWDevice


class HWDevice2(HWDevice):

    def get_properties(self):
        meta = self.init_config_info()
        for key, prop in inspect.getmembers(self.__class__):
            if isinstance(prop, property):
                meta[key] = getattr(self, key)
        return meta

    def init_config_info(self):
        return dict(name=self.name)

    def end_config_info(self, meta):
        if self.children:
            meta['children'] = [child.name for child in self.children]
        return meta


class RoboticPositioner(HWDevice2):
    def __init__(self, name, id, pos=None, parent=None):
        super(RoboticPositioner,self). __init__(name, parent)
        self._id = id
        if pos is None:
            pos = (0.0, 0.0, 0.0)
        self.x_fix = pos[0]
        self.y_fix = pos[1]
        self.pa_fix = pos[2]

        self.x_ = self.x_fix
        self.y_ = self.y_fix
        self.pa_ = self.pa_fix

        self.fb = None

    def move_to(self, x, y, pa):
        # FIME: check this is inside patrol area
        self.x_ = x
        self.y_ = y
        self.pa_ = pa

    def park(self):
        self.move_to(self.x_fix, self.y_fix, self.pa_fix)

    def connect_bundle(self, bundle):
        self.fb = bundle

    def fibers_in_focal_plane(self):
        # Given PA and center, return position
        # of the fibers
        import math
        import numpy as np
        base = math.pi / 3.0
        ini = base / 2.0
        rad = 0.5365 # Geometry
        PA = self.pa * (math.pi / 180)

        if self.fb is None:
            return ([], [])
        else:
            # To the left
            #angs = PA + ini + base * np.arange(6)
            # To the right
            angs = -PA + ini + base * np.arange(6)
            m0 = rad * np.cos(angs)
            m1 = rad * np.sin(angs)

            # Rearrange
            xx = m0[[2, 3, 1, 0, 4, 5, 0]]
            yy = m1[[2, 3, 1, 0, 4, 5, 0]]
            xx[3] = 0.0 # Center
            yy[3] = 0.0 # Center
            xx += self.x
            yy += self.y
            return self.fb.fibs_id, list(zip(xx, yy))

    @property
    def x(self):
        return self.x_

    @property
    def y(self):
        return self.y_

    @property
    def pa(self):
        return self.pa_

    def init_config_info(self):
        meta = super(RoboticPositioner, self).init_config_info()
        meta["id"] = self._id
        if self.fb:
            meta['bundle'] = self.fb.config_info()
        return meta


class FiberMOS(HWDevice2):
    def __init__(self, name):
        super(FiberMOS, self).__init__(name, parent=None)

        self.nrobots = 0

    @property
    def nbundle(self):
        return self.nrobots

    def config_info(self):
        return visit(self)


class PseudoSlit2(HWDevice):
    def __init__(self, name, insmode):

        super(PseudoSlit2, self).__init__(name)

        # Positions of the fibers in the PS slit
        self.y_pos = {}
        self.insmode = insmode

    def connect_fibers(self, fibid, pos):
        self.y_pos = dict(zip(fibid, pos))

    def config_info(self):
        return {'name': self.name, 'insmode': self.insmode}


def visit(node, root='', meta=None):
    sep = '.'
    if meta is None:
        meta = {}

    if root != '':
        node_name = root + sep + node.name
    else:
        node_name = node.name

    meta[node_name] = node.get_properties()
    submeta = meta
    for child in node.children:
        visit(child, root=node_name, meta=submeta)
    return meta