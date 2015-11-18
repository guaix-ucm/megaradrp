#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''Tests for the calibration module.'''

from megaradrp.requirements import MasterBiasRequirement, MasterDarkRequirement
from megaradrp.requirements import MasterFiberFlatRequirement


def test_MasterBiasRequirement():
    obj = MasterBiasRequirement()
    assert isinstance(obj, MasterBiasRequirement)


def test_MasterDarkRequirement():
    obj = MasterDarkRequirement()
    assert isinstance(obj, MasterDarkRequirement)


def test_MasterFiberFlatRequirement():
    obj = MasterFiberFlatRequirement()
    assert isinstance(obj, MasterFiberFlatRequirement)