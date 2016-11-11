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

"""Bad Pixel Mask (BPM) recipe"""

import uuid
import datetime

import numina.exceptions
import astropy.io.fits as fits
from numina.array.cosmetics import ccdmask
from numina.core import Product
from numina.array import combine
import numina.core.validator as val

from megaradrp.processing.combine import basic_processing_with_combination_frames
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import  MasterBPM
import megaradrp.requirements as reqs


class BadPixelsMaskRecipe(MegaraBaseRecipe):

    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    master_bpm = Product(MasterBPM)

    def __init__(self):
        super(BadPixelsMaskRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        self.logger.info('start BPM recipe')
        N = len(rinput.obresult.frames)

        self.logger.debug('we have %d images', N)
        half = N // 2
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        self.logger.debug('we have %d images', N)
        reduced1 = basic_processing_with_combination_frames(
            rinput.obresult.frames[:half],
            flow,
            method=combine.median
        )

        self.save_intermediate_img(reduced1, 'reduced1.fits')

        reduced2 = basic_processing_with_combination_frames(
            rinput.obresult.frames[half:],
            flow,
            method=combine.median
        )

        self.save_intermediate_img(reduced2, 'reduced2.fits')

        ratio, mask, sigma = ccdmask(reduced1[0].data, reduced2[0].data, mode='full')

        hdu = fits.PrimaryHDU(mask, header=reduced1[0].header)
        hdu.header['UUID'] = uuid.uuid1().hex
        hdu.header['OBJECT'] = 'MASTER BPM'
        hdu.header['IMAGETYP'] = 'BPM'
        self.set_base_headers(hdu.header)

        hdu.header['history'] = 'BPM creation time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )

        for frame in rinput.obresult.frames:
            hdu.header['history'] = "With image {}".format(
                self.datamodel.get_imgid(frame.open())
            )

        reduced = fits.HDUList([hdu])
        self.logger.info('end BPM recipe')
        return self.create_result(master_bpm=reduced)

    def validate_input(self, recipe_input):
        "Validate input of the recipe"

        obresult = recipe_input.obresult
        # Check that the number of frames is even
        nimages = len(obresult.frames)
        if nimages % 2 != 0:
            msg = 'expected even number of frames, received {} instead'.format(nimages)
            raise numina.exceptions.ValidationError(msg)
        # Continue with additional checks
        super(BadPixelsMaskRecipe, self).validate_input(recipe_input)