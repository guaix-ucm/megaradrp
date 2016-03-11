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


from .efficiency import Efficiency

class FiberBundle(object):
    def __init__(self, name, fid, bid, static=True, transmission=None, inactive=None):
        self.name = name
        # Geometry of the fibers
        self.size = 0.31
        self.fwhm = 3.6
        self.sigma = self.fwhm / 2.3548

        self.fibs_id = fid
        self.bunds_id = bid
        self.static = static
        self.nfibers = len(self.fibs_id)
        # Include the transmission of the fibers (all fibers are equal)
        if transmission is None:
            self._transmission = Efficiency()
        else:
            self._transmission = transmission

        if inactive is None:
            self.inactive = []
        else:
            self.inactive = inactive

    def transmission(self, wl):
        return self._transmission.response(wl)

    def meta(self):
        return {'name': self.name,
                'nfibers': self.nfibers,
                'fibs_id': self.fibs_id,
                'bunds_id': self.bunds_id,
                'static': self.static,
                'inactive_fibs_id': self.inactive
                }