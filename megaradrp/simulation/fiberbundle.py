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

import math

from astropy import units as u

from .efficiency import Efficiency


class FiberBundle(object):
    def __init__(self, name, bid, static=True):
        self.name = name
        # Geometry of the fibers

        self.bunds_id = bid
        self.static = static
        self.children = []

    def add_light_fiber(self, lf):
        self.children.append(lf)

    @property
    def fibs_id(self):
        return [ch.fibid for ch in self.children]

    @property
    def inactive_fibs_id(self):
        return [ch.fibid for ch in self.children if ch.inactive]

    def config_info(self):
        return {'name': self.name,
                'nfibers': len(self.children),
                'fibs_id': self.fibs_id,
                'id': self.bunds_id,
                'static': self.static,
                'inactive_fibs_id': self.inactive_fibs_id
                }