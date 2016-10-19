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
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.types import ProcessedRSS, ProcessedFrame


class LCBImageRecipe(ImageRecipe):
    """Process LCB images."""

    final = Product(ProcessedRSS)
    reduced = Product(ProcessedFrame)
    rss = Product(ProcessedRSS)
    #target = Product(ProcessedRSS)
    # sky = Product(ProcessedRSS)

    def __init__(self):
        super(LCBImageRecipe, self).__init__()

    def run(self, rinput):

        self.logger.info('starting LCB reduction')

        reduced, rss_data = super(LCBImageRecipe,self).base_run(rinput)

        return self.create_result(
            final=rss_data,
            reduced=reduced,
            rss=rss_data
        )
