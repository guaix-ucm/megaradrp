#
# Copyright 2014-2016 Universidad Complutense de Madrid
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


from numina.core import Product
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterBias
from megaradrp.requirements import MasterBPMRequirement


class BiasRecipe(MegaraBaseRecipe):
    """Process BIAS images and create MASTER_BIAS.

    Attributes
    ----------

    master_bpm: MasterBPM (requirement)
    master_bias; MasterBias (product)
    """

    master_bpm = MasterBPMRequirement()
    master_bias = Product(MasterBias)

    def run(self, rinput):
        self.logger.info('start bias recipe')
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        errors  = False
        if not errors:
            self.logger.info('not computing errors')
        hdulist = basic_processing_with_combination(rinput, flow, method=combine.median, errors=errors)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        result = self.create_result(master_bias=hdulist)
        self.logger.info('end bias recipe')
        return result
