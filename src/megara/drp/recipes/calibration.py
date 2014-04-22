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

from astropy.io import fits

from numina.core import BaseRecipe, RecipeRequirements
from numina.core import Product, DataProductRequirement
from numina.core import define_requirements, define_result
from numina.array.combine import median as c_median
#from numina.logger import log_to_history

from megara.drp.core import RecipeResult
from megara.drp.products import MasterBias, MasterDark

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
    def run(self, rinput):
        _logger.info('starting bias reduction')
        
        cdata = []
        try:
            for frame in rinput.obresult.frames:
                cdata.append(frame.open())
 
            _logger.info('stacking %d images using median', len(cdata))
            
            data = c_median([d[0].data for d in cdata], dtype='float32')
            template_header = cdata[0][0].header
            hdu = fits.PrimaryHDU(data[0], header=template_header)
        finally:
            for hdulist in cdata:
                hdulist.close()
      
        _logger.info('bias reduction ended')

        result = BiasRecipeResult(biasframe=hdu)
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
    def run(self, rinput):

        _logger.info('starting dark reduction')

        _logger.info('dark reduction ended')

        result = DarkRecipeResult(darkframe=None)
        return result

