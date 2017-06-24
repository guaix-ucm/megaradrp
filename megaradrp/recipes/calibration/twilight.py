#
# Copyright 2015-2017 Universidad Complutense de Madrid
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

""" Twilight fiber flat Calibration Recipes for Megara"""

from __future__ import division, print_function


import numpy
from astropy.io import fits

from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

import megaradrp.requirements as reqs
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterTwilightFlat
from megaradrp.types import ProcessedRSS, ProcessedFrame
# Flat 2D
from megaradrp.processing.combine import basic_processing_with_combination
from numina.array import combine
# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR, FiberFlatCorrector

import numina.core.recipeinout as recipeio
from numina.core.metarecipes import generate_docs


@generate_docs
class RecipeInput(recipeio.RecipeInput):
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_traces = reqs.MasterTraceMapRequirement()
    # master_weights = Requirement(MasterWeights, 'Set of files with extraction weights')
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()


@generate_docs
class RecipeResult(recipeio.RecipeResult):
    reduced_image = Product(ProcessedFrame)
    reduced_rss = Product(ProcessedRSS)
    master_twilightflat = Product(MasterTwilightFlat)


@recipeio.define_input(RecipeInput)
@recipeio.define_result(RecipeResult)
class TwilightFiberFlatRecipe(MegaraBaseRecipe):
    """Process TWILIGHT_FLAT images and create MASTER_TWILIGHT_FLAT product.

    This recipe process a set of continuum flat images obtained in
    **Twilight Fiber Flat** mode and returns the master twilight flat product
    The recipe also returns the result of processing the input images up to
    slitflat correction. and the result RSS of the processing
    up to wavelength calibration.

    See Also
    --------
    megaradrp.types.MasterTwilightFlat: description of MasterTwilightFlat product

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None) and corrected from pixel-to-pixel flat
    if `master_slitflat` is not None.
    Images thus corrected are the stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_traces` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    To normalize the `master_twilight_flat`, each fiber is divided by the average
    of the 200 central columns. This RSS
    image is returned in the field `master_twilightflat` of the recipe result.

    """

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        final_image = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def run(self, rinput):
        # Basic processing

        self.logger.info('twilight fiber flat reduction started')

        img = self.process_flat2d(rinput)
        self.save_intermediate_img(img, 'reduced_image.fits')

        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.master_traces, self.datamodel)
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.master_wlcalib, self.datamodel)
        flipcor = FlipLR()
        calibrator_ff = FiberFlatCorrector(rinput.master_fiberflat.open(), self.datamodel)

        img = splitter1(img)
        flat2d = splitter1.out # Copy before extraction
        img = calibrator_aper(img)
        img = splitter2(img)
        rss_base = splitter2.out # Copy before el calibration
        self.logger.debug('Flip RSS left-rigtht, before WL calibration')
        img = flipcor.run(img)
        # Calibrate in WL
        rss_wl = calibrator_wl(img)
        # Calibrate fiber flat
        rss_wl = calibrator_ff(rss_wl)

        self.save_intermediate_img(rss_wl, 'reduced_rss.fits')
        # Measure values in final
        start = 1900
        end = 2100
        self.logger.info('doing mean between columns %d-%d', start, end)
        colapse = rss_wl[:, start:end].mean(axis=1)

        normalized = numpy.tile(colapse[:, numpy.newaxis], 4096)

        template_header = flat2d[0].header
        master_t_hdu = fits.PrimaryHDU(normalized, header=template_header)
        master_t = fits.HDUList([master_t_hdu])

        self.logger.info('twilight fiber flat reduction ended')
        result = self.create_result(reduced_image=flat2d, reduced_rss=rss_wl,
                                    master_twilightflat=master_t)

        return result
