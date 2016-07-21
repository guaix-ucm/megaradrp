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


from .wheel import Carrousel


class PseudoSlit(object):
    def __init__(self, name, insmode):

        super(PseudoSlit, self).__init__()
        self.name = name
        self.insmode = insmode
        self.fiberset = None
        self.positions = []

    def connect_fibers(self, fiberset, positions):
        self.fiberset = fiberset
        self.positions = positions

    def y_pos(self, fibsid):
        result = []
        for fibid in fibsid:
            pos = self.positions[fibid-1]
            result.append(pos)
        return result

    def config_info(self):
        return {'name': self.name,
                'insmode': self.insmode,
                'nfibers': len(self.fiberset.fibers),
                }


class PseudoSlitSelector(Carrousel):
    pass