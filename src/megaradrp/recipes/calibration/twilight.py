#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Twilight fiber flat Calibration Recipes for Megara"""

import warnings

import numpy
from astropy.io import fits

from numina.core import Result, Parameter
from numina.util.flow import SerialFlow
from numina.exceptions import ValidationError

import megaradrp.requirements as reqs
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import MasterTwilightFlat
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame

# Flat 2D
from numina.processing.combine import basic_processing_with_combination
from numina.array import combine

# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import (
    FLAT_MASK_EXTNAME,
    FLAT_MASK_ZERO_ORDER,
    FLAT_MASK_ZERO_ORDER_WEAK,
    FiberFlatCorrector,
    FlipLR,
    find_flat_mask,
)


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


def pixel_2d_check_or_empty(value):
    """Check this is a valid 2D pixel list or an empty list"""
    if value is None or len(value) == 0:
        return value
    return pixel_2d_check(value)


def select_normalize_region(data, coverage, zero_order, good, region=(1900, 2100)):
    """Columns (start, end) to normalize, the window `region` moved out of the zero order

    The window starts at `region`. If a good fiber has a pixel of the zero
    order of the VPH, or a pixel without coverage, in the window, the window
    is moved column by column in the direction in which the signal grows,
    until it has neither; if there is no such window in that direction, in
    the other one. The signal of each column is the median, over the `good`
    fibers, of their pixels with coverage and out of the zero order.

    Parameters
    ----------
    data : numpy.ndarray
        RSS, shape (nfibers, npix)
    coverage : numpy.ndarray
        Boolean array of the shape of `data`, the pixels with coverage
    zero_order : numpy.ndarray
        Boolean array of the shape of `data`, the pixels in the zero order
    good : numpy.ndarray
        Boolean array (nfibers,), the fibers used
    region : tuple of int
        (start, end) of the initial window

    Returns
    -------
    tuple of int or None
        (start, end) of the window, None if there is no window
    """
    npix = data.shape[1]
    start, end = region
    width = end - start
    if width <= 0 or width > npix or not good.any():
        return None
    start = min(max(start, 0), npix - width)
    # the columns where some good fiber has zero order or no coverage
    blocked = (zero_order | ~coverage)[good].any(axis=0)
    cblock = numpy.r_[0, numpy.cumsum(blocked)]
    # number of blocked columns of the window starting at each column
    nblocked = cblock[width:] - cblock[:-width]
    if nblocked[start] == 0:
        return start, start + width
    use = (coverage & ~zero_order)[good]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        signal = numpy.nanmedian(numpy.where(use, data[good], numpy.nan), axis=0)
        before = numpy.nanmean(signal[max(start - width, 0) : start]) if start > 0 else numpy.nan
        after = numpy.nanmean(signal[start + width : start + 2 * width]) if start + width < npix else numpy.nan
    step = -1 if numpy.isfinite(before) and not (after >= before) else 1
    for direction in (step, -step):
        free = numpy.flatnonzero(nblocked[start::direction] == 0)
        if free.size > 0:
            new_start = start + direction * int(free[0])
            return new_start, new_start + width
    return None


def masked_mean(data, usable, region, good):
    """Mean of each fiber in the columns `region`, leaving out the pixels not `usable`

    The shape of the spectrum in the region is the median, over the `good`
    fibers, of each fiber divided by its mean of the usable pixels. Each fiber
    is divided by that shape, the ratio is averaged over its usable pixels and
    multiplied by the mean of the shape. Without masked pixels this is the mean
    of the fiber in the region, and a fiber with masked pixels is not biased by
    the shape of the spectrum. A fiber without usable pixels in the region gets
    the mean of all its pixels.
    """
    start, end = region
    block = data[:, start:end]
    use = usable[:, start:end]
    count = use.sum(axis=1)
    plain = block.mean(axis=1)
    rough = numpy.where(count > 0, numpy.where(use, block, 0.0).sum(axis=1) / numpy.maximum(count, 1), plain)
    norm = good & (count > 0) & (rough != 0)
    if not norm.any():
        return rough
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        shape = numpy.nanmedian(numpy.where(use[norm], block[norm] / rough[norm, numpy.newaxis], numpy.nan), axis=0)
    use = use & numpy.isfinite(shape) & (shape != 0)
    count = use.sum(axis=1)
    ratio = numpy.where(use, block / numpy.where(use, shape, 1.0), 0.0)
    return numpy.where(count > 0, ratio.sum(axis=1) / numpy.maximum(count, 1) * numpy.nanmean(shape), plain)


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

    To normalize the `master_twilight_flat`, each fiber is divided by its mean
    in the columns `normalize_region`, divided by the mean of that value over the
    good fibers (those with a wavelength solution, and not invalid in the fiber
    flat). With `normalize_region` empty (the default), the region starts at
    the columns 1900-2100, in the centre of the spectrum. If a good fiber has a
    pixel of the zero order of the VPH (marked in FLATMASK, the mask of the
    fiber flat, and not corrected in the twilight) or a pixel without coverage
    there, the region is moved in the direction in which the signal grows until
    it has neither. The pixels marked in FLATMASK (the zero order and the
    absorption windows where the fiber flat is interpolated) are left out of
    the means of `normalize_region` and `continuum_region`. This RSS image, with
    a constant value per fiber, is returned in the field `master_twilightflat` of
    the recipe result.

    """

    # the initial window of the normalization, in the centre of the spectrum
    default_normalize_region = (1900, 2100)

    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias="master_traces")
    extraction_offset = Parameter([0.0], "Offset traces for extraction", accept_scalar=True)
    normalize_region = Parameter(
        [],
        "Columns used to normalize the flat-field; [] to start at 1900-2100 and avoid the zero order",
        validator=pixel_2d_check_or_empty,
    )
    continuum_region = Parameter(
        [1900, 1900], "Subtract this region before normalize the flat-field", validator=pixel_2d_check_or_none
    )
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()

    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_twilightflat = Result(MasterTwilightFlat)

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        final_image = basic_processing_with_combination(rinput, flow, method=self.combine_median_scaled)
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

        reduced_rss = flow_1d(img)
        return reduced_rss

    def run(self, rinput):
        # Basic processing

        self.logger.info("twilight fiber flat reduction started")

        img = self.process_flat2d(rinput)
        # Copy image
        reduced_image = fits.HDUList([hdu.copy() for hdu in img])
        self.save_intermediate_img(reduced_image, "reduced_image.fits")

        reduced_rss = self.run_reduction_1d(
            img,
            rinput.master_apertures,
            rinput.master_wlcalib,
            rinput.master_fiberflat,
            offset=rinput.extraction_offset,
        )

        self.save_intermediate_img(reduced_rss, "reduced_rss.fits")
        # Measure values in final
        rss_wl_data = reduced_rss[0].data
        mask = self.good_ids_mask(rinput.master_wlcalib)
        invalid = self.invalid_fibers(reduced_rss)
        if invalid:
            self.logger.info("fibers invalid in the fiber flat, not used in the normalization: %s", invalid)
            mask[[fibid - 1 for fibid in invalid if fibid <= mask.size]] = False

        shape = rss_wl_data.shape
        coverage = reduced_rss["WLMAP"].data > 0 if "WLMAP" in reduced_rss else numpy.ones(shape, dtype=bool)
        flat_mask = find_flat_mask(reduced_rss)
        flat_mask_data = flat_mask.data if flat_mask is not None else numpy.zeros(shape, dtype="uint8")
        masked = flat_mask_data > 0
        zero_order = (flat_mask_data & (FLAT_MASK_ZERO_ORDER | FLAT_MASK_ZERO_ORDER_WEAK)) > 0
        usable = coverage & ~masked

        if rinput.normalize_region:
            start, end = rinput.normalize_region
        else:
            start, end = self.default_normalize_region
            region = select_normalize_region(rss_wl_data, coverage, zero_order, mask, region=(start, end))
            if region is None:
                self.logger.warning(
                    "no region without zero order and with coverage in all the fibers, using %d-%d", start, end
                )
            else:
                if region != (start, end):
                    self.logger.info(
                        "columns %d-%d with zero order or without coverage, moved to %d-%d", start, end, *region
                    )
                start, end = region
        self.logger.info(
            "doing mean between columns %d-%d, %d pixels masked in the good fibers",
            start,
            end,
            numpy.count_nonzero(masked[mask, start:end]),
        )
        collapse = masked_mean(rss_wl_data, usable, (start, end), mask)
        if rinput.continuum_region is not None:
            start_c, end_c = rinput.continuum_region
            if end_c > start_c:
                self.logger.info("subtract mean of columns %d-%d", start_c, end_c)
                collapse_c = masked_mean(rss_wl_data, usable, (start_c, end_c), mask)
                collapse -= collapse_c

        # Normalize the collapsed array
        collapse_good = collapse[mask]
        collapse_norm = collapse / collapse_good.mean()
        normalized = numpy.tile(collapse_norm[:, numpy.newaxis], rss_wl_data.shape[1])

        # the mask of the fiber flat of the reduced RSS does not apply to the constant values of the product
        master_t = fits.HDUList([hdu.copy() for hdu in reduced_rss if hdu.name != FLAT_MASK_EXTNAME])
        master_t[0].data = normalized
        self.set_product_headers(master_t[0].header)
        master_t[0].header["history"] = f"Normalized with the mean of the columns {start}-{end}"

        self.logger.info("twilight fiber flat reduction ended")
        result = self.create_result(reduced_image=reduced_image, reduced_rss=reduced_rss, master_twilightflat=master_t)

        return result

    @staticmethod
    def invalid_fibers(rss):
        """The fibers marked invalid in the FIBERS extension of `rss`"""
        if "FIBERS" not in rss:
            return []
        hdr = rss["FIBERS"].header
        return sorted(
            int(key[3:6])
            for key, value in hdr.items()
            if key.startswith("FIB") and key.endswith("_V") and value is False
        )

    def good_ids_mask(self, calibration):

        # Bad fibers, join in a new list, the calibration is not modified
        bad_fibers = list(calibration.missing_fibers) + list(calibration.error_fitting)

        bad_idxs = [fibid - 1 for fibid in bad_fibers]

        good_idxs_mask = numpy.ones((calibration.total_fibers,), dtype="bool")
        good_idxs_mask[bad_idxs] = False
        return good_idxs_mask

    def set_product_headers(self, hdr):
        """Set metadata in the FITS headers of the main product."""
        hdr = self.set_base_headers(hdr)
        hdr["NUMTYPE"] = ("MasterTwilightFlat", "Product type")
        return hdr

    def combine_median_scaled(self, arrays, masks=None, dtype=None, out=None, zeros=None, scales=None, weights=None):

        median_vals = numpy.array([numpy.median(arr) for arr in arrays])
        self.logger.info("median values are %s", median_vals)
        # normalize by max value
        median_max = numpy.max(median_vals)
        scales = median_vals / median_max
        self.logger.info("scale values are %s", scales)
        return combine.median(arrays, masks=masks, scales=scales, dtype=dtype)


# The input and result classes, with their names before the fields were
# declared in the recipe
RecipeInput = TwilightFiberFlatRecipe.RecipeInput
RecipeResult = TwilightFiberFlatRecipe.RecipeResult
