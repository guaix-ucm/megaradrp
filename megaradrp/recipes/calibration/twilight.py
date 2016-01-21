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

""" Twilight fiber flat Calibration Recipes for Megara"""

from __future__ import division, print_function

from numina.core import Product
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterFiberFlat
from megaradrp.requirements import MasterBiasRequirement


class TwilightFiberFlatRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(TwilightFiberFlatRecipe, self).__init__(
            version="0.1.0"
        )

    def run(self, rinput):
        pass
