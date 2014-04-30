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

from numina.core import RecipeRequirements
from numina.core.recipes import BaseRecipe
from numina.core import Requirement, Product, DataProductRequirement
from numina.core.requirements import ObservationResultRequirement
from numina.core import define_requirements, define_result

from megara.drp.core import RecipeResult
from megara.drp.products import MasterBias

_logger = logging.getLogger('numina.recipes.megara')

class BiasRecipeRequirements(RecipeRequirements):
    obs = ObservationResultRequirement()

class BiasRecipeResult(RecipeResult):
    biasframe = Product(MasterBias)

@define_requirements(BiasRecipeRequirements)
@define_result(BiasRecipeResult)
class BiasRecipe(BaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    def __init__(self):
        super(BiasRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, recipe_input):
        _logger.info('starting bias reduction')

        _logger.info('stacking images')
        _logger.info('bias reduction ended')

        result = BiasRecipeResult(biasframe=None)
        return result

