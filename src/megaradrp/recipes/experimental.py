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

'''Recipes with experimental features.'''

import logging


from numina.core.recipes import BaseRecipeAutoQC
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

from megaradrp.products import MasterBias

_logger = logging.getLogger('numina.recipes.megara')


class BiasRecipe(BaseRecipeAutoQC):
    '''Process BIAS images and create MASTER_BIAS.'''

    obs = ObservationResultRequirement()

    biasframe = Product(MasterBias)

    def __init__(self):
        super(BiasRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, recipe_input):
        _logger.info('starting bias reduction')

        _logger.info('stacking images')
        _logger.info('bias reduction ended')

        result = self.create_result(biasframe=None)
        return result
