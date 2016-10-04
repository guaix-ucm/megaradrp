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

"""Calibration Recipes for Megara"""

import numpy
from numina.core import Product

from megaradrp.processing.datamodel import MegaraDataModel
from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import MasterFiberFlat
from megaradrp.processing.wavecalibration import WavelengthCalibrator


class LCBImageRecipe(ImageRecipe):
    """Process LCB images."""

    final = Product(MasterFiberFlat)
    reduced = Product(MasterFiberFlat)
    rss = Product(MasterFiberFlat)
    #target = Product(MasterFiberFlat)
    # sky = Product(MasterFiberFlat)

    def __init__(self):
        super(LCBImageRecipe, self).__init__()

    def run(self, rinput):

        self.logger.info('starting LCB reduction')

        reduced, rss_data = super(LCBImageRecipe,self).run(rinput)

        # FIXME: Flip L-R image before calibrating WL
        # Eventually this should not be necessary

        self.logger.debug('Flip RSS left-rigtht, before WL calibration')
        rss_data[0].data = numpy.fliplr(rss_data[0].data)

        datamodel = MegaraDataModel()
        calibrator = WavelengthCalibrator(rinput.wlcalib, datamodel)

        rss_wl = calibrator(rss_data)

        return self.create_result(
            final=rss_wl,
            reduced=reduced,
            rss=rss_data
        )
        # return self.create_result(final=rss, target=reduced, sky=reduced)
