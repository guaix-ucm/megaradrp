#
# Copyright 2015 Universidad Complutense de Madrid
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

"""Typical requirements of recipes"""

from numina.core import Requirement

from .products import MasterBias, MasterDark, MasterFiberFlat, MasterBPM


class MasterBiasRequirement(Requirement):
    def __init__(self, optional=False):
        super(MasterBiasRequirement,
              self).__init__(MasterBias,
                             'Master BIAS image',
                             optional=optional
                             )


class MasterBPMRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBPMRequirement,
              self).__init__(MasterBPM,
                             'Master Bad Pixel Mask',
                             optional=optional
                             )


class MasterDarkRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterDarkRequirement,
              self).__init__(MasterDark, 'Master DARK image',
                             optional=optional)


class MasterFiberFlatRequirement(Requirement):
    def __init__(self):
        super(MasterFiberFlatRequirement,
              self).__init__(MasterFiberFlat,
                             'Master fiber flat calibration'
                             )
