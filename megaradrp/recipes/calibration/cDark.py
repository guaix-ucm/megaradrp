__author__ = 'Pica4x6'

import logging

from numina.core import Product

from cBase import BaseRecipeAutoQC as MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.products import MasterDark

_logger = logging.getLogger('numina.recipes.megara')

class DarkRecipe(MegaraBaseRecipe):

    '''Process DARK images and provide MASTER_DARK. '''

    # master_bias = MasterBiasRequirement()

    darkframe = Product(MasterDark)

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

        result = self.create_result(darkframe=None)
        return result