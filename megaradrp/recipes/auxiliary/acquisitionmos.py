#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Acquire with MOS Recipe for Megara"""


from __future__ import division, print_function

import logging


from numina.core import RecipeError

from megaradrp.core.recipe import MegaraBaseRecipe


#from astropy.io import fits

from numina.core import Product, DataFrameType
from numina.core.requirements import ObservationResultRequirement, Requirement
from numina.array.combine import median

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterFiberFlat, WavelengthCalibration
from megaradrp.products import MasterWeights
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement, MasterFiberFlatRequirement
from megaradrp.requirements import MasterSlitFlatRequirement, MasterTwilightRequirement
from megaradrp.requirements import MasterTraceMapRequirement
from megaradrp.processing.fiberflat import FiberFlatCorrector
from megaradrp.processing.twilight import TwilightCorrector
from megaradrp.processing.weights import WeightsCorrector
from megaradrp.processing.combine import basic_processing_with_combination

_logger = logging.getLogger('numina.recipes.megara')


class AcquireMOSRecipe(MegaraBaseRecipe):
    """Process Focus images and find best focus."""

    obresult = ObservationResultRequirement()
    master_bpm = MasterBPMRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_slitflat = MasterSlitFlatRequirement(optional=True)

    #master_tracemap = MasterTraceMapRequirement()
    #master_fiberflat = MasterFiberFlatRequirement()
    #master_wlcalib = Requirement(WavelengthCalibration, 'Wavelength calibration table')
    # master_weights = Requirement(MasterWeights, 'Set of files')

    #master_twilight = MasterTwilightRequirement()

    frame = Product(DataFrameType)

    def __init__(self):
        super(AcquireMOSRecipe, self).__init__("1")

    def run(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration.values)

        hdulist = basic_processing_with_combination(rinput, flow, method=median)

        return self.create_result(frame=hdulist)
