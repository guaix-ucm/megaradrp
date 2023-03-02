#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

import numpy
from numina.instrument.hwdevice import HWDevice

from megaradrp.instrument.focalplane import TargetType


class RoboticPositioner(HWDevice):

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

        self._target_priority = 0
        self._target_type = TargetType.UNASSIGNED
        self._target_name = 'unknown'

        self.fb = None
        self.rpatrol = 10.0

    def move_to(self, x, y, pa):
        dis  = math.hypot(x - self.x_fix, y - self.y_fix)
        if dis >= self.rpatrol:
            # impossible to move
            raise ValueError('movement out of patrol area')
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
            angs = -PA + ini + base * numpy.arange(6)
            m0 = rad * numpy.cos(angs)
            m1 = rad * numpy.sin(angs)

            # Rearrange, fibers are counted in different order
            xx = m0[[2, 3, 1, 0, 4, 5, 0]]
            yy = m1[[2, 3, 1, 0, 4, 5, 0]]
            # This is the central fiber
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

    @property
    def position(self):
        return self.x_, self.y_, self.pa_

    @position.setter
    def position(self, value):
        self.move_to(value[0], value[1], value[2])

    @property
    def position_relative(self):
        return self.x_ - self.x_fix, self.y_ - self.y_fix, self.pa_

    @position_relative.setter
    def position_relative(self, value):
        self.move_to(value[0] + self.x_, value[1] + self.y_, value[2] + self.pa_)

    @property
    def target_priority(self):
        return self._target_priority

    @target_priority.setter
    def target_priority(self, value):
        self._target_priority = value

    @property
    def target_type(self):
        return self._target_type

    @target_type.setter
    def target_type(self, value):
        self._target_type = value

    @property
    def target_name(self):
        return self._target_name

    @target_name.setter
    def target_name(self, value):
        self._target_name = value

    def init_config_info(self):
        meta = super(RoboticPositioner, self).init_config_info()
        meta["id"] = self._id
        if self.fb:
            meta['bundle'] = self.fb.config_info()
            _, pos = self.fibers_in_focal_plane()
            meta['fibers_pos'] = pos
        return meta


class BaseFibersPlane(HWDevice):
    def __init__(self, name, fiberset, parent=None):

        super(BaseFibersPlane, self).__init__(name, parent)
        self.fiberset = fiberset
        self.confid_ = 1

    @property
    def nbundles(self):
        return self.nfibers // 7

    @property
    def nfibers(self):
        return len(self.fiberset.fibers)

    @property
    def size(self):
        return self.fiberset.size

    @property
    def sigma(self):
        return self.fiberset.sigma

    @property
    def area(self):
        return self.fiberset.area

    @property
    def conf_id(self):
        return self.confid_

    @conf_id.setter
    def conf_id(self, value):
        self.confid_ = value

    def transmission(self, wlin):
        return self.fiberset.transmission(wlin)


class FiberMOS(BaseFibersPlane):
    def __init__(self, name, fiberset, parent=None):
        super(FiberMOS, self).__init__(name, fiberset, parent)

    def fibers_in_focal_plane(self):
        fibid = []
        pos = []
        for robot in self.children.values():
            # Compute positions of the fibers
            # in focal plane for fibers
            res1, res2 = robot.fibers_in_focal_plane()
            fibid.extend(res1)
            pos.extend(res2)
        # This must be sorted
        fibid = numpy.array(fibid)
        pos = numpy.array(pos)
        idx = numpy.argsort(fibid)
        return fibid[idx], pos[idx]

    def robots_in_focal_plane(self):
        robid = []
        pos = []
        for robot in self.children.values():
            # Compute positions of the fibers
            # in focal plane for fibers
            res2 = robot.position
            robid.append(robot._id)
            pos.append(res2[:2])
        return robid, pos


class LargeCompactBundle(BaseFibersPlane):
    def __init__(self, name, fiberset, lcb_pos, parent=None):
        super(LargeCompactBundle, self).__init__(name, fiberset, parent)
        self.lcb_pos = lcb_pos

    def fibers_in_focal_plane(self):
        fibid = [fiber.fibid for fiber in self.fiberset.fibers.values()]
        pos = [self.lcb_pos[fiber.fibid] for fiber in self.fiberset.fibers.values()]

        return fibid, pos

    @property
    def fibers(self):
        res = []
        for bundid in self.fiberset.bundles:
            for lf in self.fiberset.bundles[bundid].lf:
                pos = self.lcb_pos[lf.fibid]
                res.append((lf.fibid, bundid, pos[0], pos[1], lf.inactive))

        return res
