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

from megaradrp.types import ProcessedMultiRSS
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.processing.multirss import generate_multi_rss


class LCBFastMappingRecipe(MegaraBaseRecipe):
    """Process LCB Fast Mapping Recipe."""

    final = Product(ProcessedMultiRSS)

    def run(self, rinput):
        self.logger.info('start FastMappingRecipe')
        obresult = rinput.obresult
        imgs = [frame.open() for frame in obresult.frames]

        result = generate_multi_rss(imgs)
        self.logger.info('end FastMappingRecipe')

        return self.create_result(final=result)
