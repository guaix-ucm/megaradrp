from __future__ import division, print_function

import logging

from megaradrp.products import MasterFiberFlat
from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement

from numina.core import Product
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement


_logger = logging.getLogger('numina.recipes.megara')


class TwilightFiberFlatRecipe(MegaraBaseRecipe):

    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()

    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(ArrayType)

    def __init__(self):
        super(TwilightFiberFlatRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        pass