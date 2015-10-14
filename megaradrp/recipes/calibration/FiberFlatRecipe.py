from __future__ import division, print_function

import logging

from astropy.io import fits

from megaradrp.core import apextract_tracemap
from megaradrp.products import TraceMap, MasterFiberFlat
from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement

from numina.core import Requirement, Product
from numina.core.requirements import ObservationResultRequirement

_logger = logging.getLogger('numina.recipes.megara')


class FiberFlatRecipe(MegaraBaseRecipe):
    '''Process FIBER_FLAT images and create MASTER_FIBER_FLAT.'''

    # Requirements
    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
    # Products
    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        # Basic processing
        reduced = self.bias_process_common(rinput.obresult, rinput.master_bias)

        _logger.info('extract fibers')
        rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        # FIXME: we are ignoring here all the possible bad pixels
        # and WL distortion when doing the normalization
        # rssdata /= rssdata.mean() #Originally uncomment
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        rss = fits.HDUList([rsshdu])

        _logger.info('extraction completed')
        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced, fiberflat_rss=rss)
        return result