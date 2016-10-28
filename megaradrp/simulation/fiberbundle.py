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


class FiberBundle(object):
    def __init__(self, name, bid, static=True):
        super(FiberBundle, self).__init__()
        # Geometry of the fibers
        self.name = name
        self.bunds_id = bid
        self.static = static
        self.lf = []

    def add_light_fiber(self, lf):
        self.lf.append(lf)

    @property
    def fibs_id(self):
        return [ch.fibid for ch in self.lf]

    @property
    def inactive_fibs_id(self):
        return [ch.fibid for ch in self.lf if ch.inactive]

    @property
    def active_fibs_id(self):
        return [ch.fibid for ch in self.lf if not ch.inactive]

    def config_info(self):
        return {'name': self.name,
                'nfibers': len(self.lf),
                'fibs_id': self.fibs_id,
                'id': self.bunds_id,
                'static': self.static,
                'inactive_fibs_id': self.inactive_fibs_id
                }