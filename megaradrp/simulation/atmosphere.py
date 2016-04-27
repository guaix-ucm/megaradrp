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

from megaradrp.simulation.efficiency import InterpolFile

class AtmosphereModel(object):

    def __init__(self, twilight, nightsky):
        self.tw_interp = twilight
        self.ng_interp = nightsky

    def twilight_spectrum(self, wl_in):
        """Twilight spectrum"""
        return 5e4 * self.tw_interp(wl_in)

    def night_spectrum(self, wl_in):
        """Night spectrum"""
        return 5e3 * self.ng_interp(wl_in)
