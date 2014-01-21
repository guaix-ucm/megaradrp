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

'''Calibration Recipes for Megara'''

import logging

from numina.core import BaseRecipe, RecipeRequirements
from numina.core import Requirement, Product, DataProductRequirement
from numina.core import define_requirements, define_result
#from numina.logger import log_to_history

from megara.core import RecipeResult
from megara.products import MasterBias, MasterDark

_logger = logging.getLogger('numina.recipes.megara')

class BiasRecipeRequirements(RecipeRequirements):
    pass

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

    # FIXME find a better way of doing this automatically
    # @log_to_history(_logger)
    def run(self, obresult, reqs):
        _logger.info('starting bias reduction')

        _logger.info('stacking images')
        _logger.info('bias reduction ended')

        result = BiasRecipeResult(biasframe=None)
        return result

class DarkRecipeRequirements(BiasRecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')

class DarkRecipeResult(RecipeResult):
    darkframe = Product(MasterDark)

@define_requirements(DarkRecipeRequirements)
@define_result(DarkRecipeResult)
class DarkRecipe(BaseRecipe):
    '''Process DARK images and provide MASTER_DARK. '''

    def __init__(self):
        super(DarkRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    # FIXME find a better way of doing this automatically
    # @log_to_history(_logger)
    def run(self, obresult, reqs):

        _logger.info('starting dark reduction')

        _logger.info('dark reduction ended')

        result = DarkRecipeResult(darkframe=None)
        return result

