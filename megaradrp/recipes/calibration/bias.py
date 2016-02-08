#
# Copyright 2014-2015 Universidad Complutense de Madrid
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

import logging

from astropy.io import fits

from numina.core import Product, RecipeError
from numina.core.requirements import ObservationResultRequirement

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterBias
from megaradrp.requirements import MasterBPMRequirement

_logger = logging.getLogger('numina.recipes.megara')


class BiasRecipe(MegaraBaseRecipe):
    """Process BIAS images and create MASTER_BIAS."""
    obresult = ObservationResultRequirement()

    master_bias = Product(MasterBias)

    def __init__(self):
        super(BiasRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):

        _logger.info('starting bias reduction')

        if not rinput.obresult.images:
            raise RecipeError('Frame list is empty')

        hdu, data = self.hdu_creation(rinput.obresult)

        hdr = hdu[0].header
        hdr = self.set_base_headers(hdr)
        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        hdulist = fits.HDUList(hdu + [varhdu, num])
        _logger.info('bias reduction ended')

        result = self.create_result(master_bias=hdulist)
        return result
