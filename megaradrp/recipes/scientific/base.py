#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

from numina.core import Product
from numina.core.requirements import ObservationResultRequirement, Requirement
# from numina.flow import SerialFlow

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterFiberFlat, WavelengthCalibration
from megaradrp.products import MasterWeights, TraceMap
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement, MasterFiberFlatRequirement
from megaradrp.requirements import MasterSlitFlatRequirement, MasterTwilightRequirement
# from megaradrp.processing.fiberflat import FiberFlatCorrector
# from megaradrp.processing.twilight import TwilightCorrector
# from megaradrp.processing.weights import WeightsCorrector
from megaradrp.core.processing import apextract_tracemap_2

_logger = logging.getLogger('numina.recipes.megara')


class ImageRecipe(MegaraBaseRecipe):
    """Base Image."""

    # Requirements  
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    master_slitflat = MasterSlitFlatRequirement()
    wlcalib = Requirement(WavelengthCalibration, 'Wavelength calibration table')
    # master_weights = Requirement(MasterWeights, 'Set of files')
    master_fiberflat = MasterFiberFlatRequirement()
    master_twilight = MasterTwilightRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')

    # Products
    # final = Product(MasterFiberFlat)
    # target = Product(MasterFiberFlat)
    # sky = Product(MasterFiberFlat)

    def __init__(self):
        super(ImageRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):

        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)
        rssdata = apextract_tracemap_2(reduced[0].data, rinput.tracemap)
        # _logger.info('Starting: resample_rss_flux')
        # final, wcsdata = self.resample_rss_flux(reduced[0].data, self.get_wlcalib(rinput.wlcalib))
        #
        # # Add WCS spectral keywords
        # hdu_f = fits.PrimaryHDU(final, header=reduced[0].header)
        # master_lcb = fits.HDUList([hdu_f])
        #
        # return master_lcb, reduced, reduced
        return reduced, rssdata