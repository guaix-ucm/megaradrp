#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

""" Twilight fiber flat Calibration Recipes for Megara"""

from __future__ import division, print_function


import numpy
from astropy.io import fits

from numina.core import Result, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.util.flow import SerialFlow
from numina.exceptions import ValidationError

import megaradrp.requirements as reqs
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import MasterTwilightFlat
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame
# Flat 2D
from megaradrp.processing.combine import basic_processing_with_combination
from numina.array import combine
# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR, FiberFlatCorrector

import numina.core.recipeinout as recipeio
from numina.core.metarecipes import generate_docs


def pixel_2d_check(value):
    """Check this is a valid 2D pixel list"""
    if len(value) != 2:
        raise ValidationError("must have 2 elements")

    return value


def pixel_2d_check_or_none(value):
    """Check this is a valid 2D pixel list or None"""
    if value is None:
        return value
    return pixel_2d_check(value)


@generate_docs
class RecipeInput(recipeio.RecipeInput):
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias='master_traces')
    extraction_offset = Parameter([0.0], 'Offset traces for extraction', accept_scalar=True)
    normalize_region = Parameter([1900, 2100], 'Region used to normalize the flat-field',
                                 validator=pixel_2d_check)
    continuum_region = Parameter([1900, 1900], 'Subtract this region before normalize the flat-field',
                                 validator=pixel_2d_check_or_none)
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()


@generate_docs
class RecipeResult(recipeio.RecipeResult):
    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_twilightflat = Result(MasterTwilightFlat)


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
    `master_apertures` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    To normalize the `master_twilight_flat`, each fiber is divided by the average
    of the column range given in `normalize_region`. This RSS
    image is returned in the field `master_twilightflat` of the recipe result.

    """

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        final_image = basic_processing_with_combination(
            rinput, flow, method=self.combine_median_scaled
        )
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def run_reduction_1d(self, img, tracemap, wlcalib, fiberflat, offset=None):
        # 1D, extraction, Wl calibration, Flat fielding
        correctors = []
        correctors.append(ApertureExtractor(tracemap, self.datamodel, offset=offset))
        correctors.append(FlipLR())
        correctors.append(WavelengthCalibrator(wlcalib, self.datamodel))
        correctors.append(FiberFlatCorrector(fiberflat.open(), self.datamodel))

        flow_1d = SerialFlow(correctors)

        reduced_rss =  flow_1d(img)
        return reduced_rss

    def run(self, rinput):
        # Basic processing

        self.logger.info('twilight fiber flat reduction started')

        img = self.process_flat2d(rinput)
        # Copy image
        reduced_image = fits.HDUList([hdu.copy() for hdu in img])
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        reduced_rss = self.run_reduction_1d(img,
                                            rinput.master_apertures,
                                            rinput.master_wlcalib,
                                            rinput.master_fiberflat,
                                            offset=rinput.extraction_offset
                                            )

        self.save_intermediate_img(reduced_rss, 'reduced_rss.fits')
        # Measure values in final
        rss_wl_data = reduced_rss[0].data
        mask = self.good_ids_mask(rinput.master_wlcalib)

        start, end = rinput.normalize_region
        self.logger.info('doing mean between columns %d-%d', start, end)
        colapse = rss_wl_data[:, start:end].mean(axis=1)
        if rinput.continuum_region is not None:
            start_c, end_c = rinput.continuum_region
            if end_c > start_c:
                self.logger.info('subtract mean of columns %d-%d', start_c, end_c)
                colapse_c = rss_wl_data[:, start_c:end_c].mean(axis=1)
                colapse -= colapse_c

        # Normalize the colapsed array
        colapse_good = colapse[mask]
        colapse_norm = colapse / colapse_good.mean()
        normalized = numpy.tile(colapse_norm[:, numpy.newaxis], rss_wl_data.shape[1])

        master_t = fits.HDUList([hdu.copy() for hdu in reduced_rss])
        master_t[0].data = normalized
        self.set_base_headers(master_t[0].header)

        self.logger.info('twilight fiber flat reduction ended')
        result = self.create_result(reduced_image=reduced_image,
                                    reduced_rss=reduced_rss,
                                    master_twilightflat=master_t
                                    )

        return result

    def good_ids_mask(self, calibration):

        # Bad fibers, join:
        bad_fibers = calibration.missing_fibers
        bad_fibers.extend(calibration.error_fitting)

        bad_idxs = [fibid - 1 for fibid in bad_fibers]

        good_idxs_mask = numpy.ones((calibration.total_fibers,), dtype='bool')
        good_idxs_mask[bad_idxs] = False
        return good_idxs_mask

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(TwilightFiberFlatRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterTwilightFlat', 'Product type')
        return hdr

    def combine_median_scaled(self, arrays, masks=None, dtype=None, out=None,
                                    zeros=None, scales=None,
                                    weights=None):

        median_vals = numpy.array([numpy.median(arr) for arr in arrays])
        self.logger.info("median values are %s", median_vals)
        # normalize by max value
        median_max = numpy.max(median_vals)
        scales = median_max / median_vals
        self.logger.info("scale values are %s", scales)
        return combine.median(arrays, scales=scales, dtype=dtype)
