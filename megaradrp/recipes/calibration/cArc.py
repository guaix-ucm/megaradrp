__author__ = 'Pica4x6'

import logging

from numina.core import Product
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement

from cBase import BaseRecipeAutoQC as MegaraBaseRecipe

from megaradrp.products import MasterFiberFlat
from megaradrp.requirements import MasterBiasRequirement

_logger = logging.getLogger('numina.recipes.megara')

class ArcRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(ArcRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        pass