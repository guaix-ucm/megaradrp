#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

"""Fiber flat calibration Recipe for Megara"""

from __future__ import division, print_function

import logging

from astropy.io import fits

from numina.core import Product, Requirement
from numina.core.requirements import ObservationResultRequirement

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterFiberFlat
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement
from megaradrp.core.processing import apextract_tracemap

_logger = logging.getLogger('numina.recipes.megara')


class FiberFlatRecipe(MegaraBaseRecipe):
    """Process FIBER_FLAT images and create MASTER_FIBER_FLAT."""

    # Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')
    # Products
    fiberflat_frame = Product(MasterFiberFlat)
    master_fiberflat = Product(MasterFiberFlat)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        # Basic processing
        parameters = self.get_parameters(rinput)

        reduced = self.bias_process_common(rinput.obresult, parameters)

        _logger.info('extract fibers')
        rssdata = apextract_tracemap(reduced[0].data, rinput.tracemap)
        # FIXME: we are ignoring here all the possible bad pixels
        # and WL distortion when doing the normalization
        # rssdata /= rssdata.mean() #Originally uncomment
        rsshdu = fits.PrimaryHDU(rssdata, header=reduced[0].header)
        rss = fits.HDUList([rsshdu])

        _logger.info('extraction completed')
        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced, master_fiberflat=rss)
        return result