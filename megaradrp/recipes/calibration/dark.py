#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

"""Dark current calibration Recipe for Megara"""

from numina.core import Product
from numina.array import combine

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.types import MasterDark
from megaradrp.processing.combine import basic_processing_with_combination


class DarkRecipe(MegaraBaseRecipe):

    """Process DARK images and provide MASTER_DARK."""

    master_bias = MasterBiasRequirement()

    master_dark = Product(MasterDark)

    def __init__(self):
        super(DarkRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):

        flow = self.init_filters(rinput, rinput.obresult.configuration.values)
        hdulist = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        result = self.create_result(master_dark=hdulist)
        return result
