#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Bad Pixel Mask (BPM) recipe"""

import uuid
import datetime

import numina.exceptions
import astropy.io.fits as fits
from numina.array.cosmetics import ccdmask
from numina.core import Result
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination_frames
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import  MasterBPM
import megaradrp.requirements as reqs


class BadPixelsMaskRecipe(MegaraBaseRecipe):
    """Process defocussed FIBER_FLAT images and create MASTER_BPM product.

    This recipe process a set of defocused continuum flat images obtained in
    **Bad-pixels mask** mode and returns the master bad-pixels mask product.

    See Also
    --------
    numina.array.cosmetics.ccdmask: algorithm to select bad-pixels
    megaradrp.types.MasterBPM: description of MasterBPM product

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bias and dark current (if `master_dark` is not None). The first
    half of the images are the stacked using the median and saved as
    intermediate result 'reduced_image_1.fits'. The second half is also
    combined and saved as intermediate result 'reduced_image_2.fits'

    These two images are passed to the `ccdmask` function, that selects
    bad-pixels by finding outliers in the ratio of the two images.

    The mask is returned in the field `master_bpm` of the recipe result.
    """

    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    master_bpm = Result(MasterBPM)

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

        self.save_intermediate_img(reduced1, 'reduced_image_1.fits')

        reduced2 = basic_processing_with_combination_frames(
            rinput.obresult.frames[half:],
            flow,
            method=combine.median
        )

        self.save_intermediate_img(reduced2, 'reduced_image_2.fits')

        ratio, mask, sigma = ccdmask(reduced1[0].data, reduced2[0].data, mode='full')

        hdu = fits.PrimaryHDU(mask, header=reduced1[0].header)
        hdu.header['UUID'] = str(uuid.uuid1())

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
        """Validate input of the recipe.

        The number of frames in `recipe_input.obresult` must be even.

        Raises
        ------
        numina.exceptions.ValidationError
            If the number of frames in `obresult` is odd

        """

        obresult = recipe_input.obresult
        # Check that the number of frames is even
        nimages = len(obresult.frames)
        if nimages % 2 != 0:
            msg = 'expected even number of frames, received {} instead'.format(nimages)
            raise numina.exceptions.ValidationError(msg)
        # Continue with additional checks
        return super(BadPixelsMaskRecipe, self).validate_input(recipe_input)

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(BadPixelsMaskRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterBPM', 'Product type')
        return hdr