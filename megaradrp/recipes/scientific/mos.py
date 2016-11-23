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

from numina.core import Product

from megaradrp.types import ProcessedRSS, ProcessedFrame
from .base import ImageRecipe


class MOSImageRecipe(ImageRecipe):
    """Process MOS images."""

    reduced = Product(ProcessedFrame)
    final = Product(ProcessedRSS)
    target = Product(ProcessedRSS)
    sky = Product(ProcessedRSS)

    def run(self, rinput):
        self.logger.info('starting MOS reduction')

        reduced2d, rss_data = super(MOSImageRecipe, self).base_run(rinput)

        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data)
        self.logger.info('end sky subtraction')
        self.logger.info('end MOS reduction')

        return self.create_result(
            reduced=reduced2d,
            final=final,
            target=origin,
            sky=sky
        )
