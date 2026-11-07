#
# Copyright 2011-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Fiber flat calibration Recipe for Megara"""

import numpy
from astropy.io import fits
from numina.core import Result, Parameter
import numina.exceptions

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.ntypes import MasterFiberFlat
import megaradrp.requirements as reqs
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame

# Flat 2D
from numina.processing.combine import basic_processing_with_combination
from numina.array import combine

# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR


def _smoothing_window_check(val):
    """Check 'smoothing_window'"""
    # raise ValueError("something")
    if val < 5:
        raise numina.exceptions.ValidationError("must be >= 5")
    if val % 2 == 0:
        raise numina.exceptions.ValidationError("must be odd")
    return val


class FiberFlatRecipe(MegaraBaseRecipe):
    """Process FIBER_FLAT images and create MASTER_FIBER_FLAT product.

    This recipe process a set of continuum flat images obtained in
    **Fiber Flat** mode and returns the master fiber flat product
    The recipe also returns the result of processing the input images up to
    slitflat correction. and the result RSS of the processing
    up to wavelength calibration.

    See Also
    --------
    megaradrp.products.MasterFiberFlat: description of MasterFiberFlat product
    megaradrp.processing.aperture: aperture extraction
    megaradrp.processing.wavecalibration: resampling for wavelength calibration

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
    `master_apertures` and resampled accoding to the wavelength calibration in
    `master_wlcalib`. The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    To normalize the `master_fiberflat`, each fiber is divided by a smoothed
    version (using a Savitzky-Golay filter) of the average of the valid fibers.
    Finally, all the pixels with information are fiiled with ones. This RSS
    image is returned in the field `master_fiberflat` of the recipe result.

    """

    # Requirements
    method = Parameter("median", description="Combination method", choices=["mean", "median", "mediancr", "sigmaclip"])
    method_kwargs = Parameter(dict(), description="Arguments for combination method", optional=True)
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias="master_traces")
    smoothing_window = Parameter(31, "Window for smoothing (must be odd)", validator=_smoothing_window_check)
    smoothing_knot_step = Parameter(
        150.0, "Distance in Angstrom between the knots of the spline that smooths the average spectrum of the fibers"
    )
    extraction_offset = Parameter([0.0], "Offset traces for extraction", accept_scalar=True)
    master_wlcalib = reqs.WavelengthCalibrationRequirement()

    # Results
    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_fiberflat = Result(MasterFiberFlat)

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        fmethod = getattr(combine, rinput.method)
        final_image = basic_processing_with_combination(
            rinput,
            flow,
            method=fmethod,
            method_kwargs=rinput.method_kwargs,
        )
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def obtain_fiber_flat(
        self, rss_wl, col1=1900, col2=2100, window=31, degree=3, knot_step=150.0, min_level=0.3, min_coverage=0.2
    ):
        """Divide each fiber by the smoothed average spectrum of the fibers

        Each valid fiber is normalized by its median in the columns `col1:col2`.
        The fibers fainter than `min_level` times the median of the fibers are
        left out. In each wavelength, the average spectrum is the mean of the
        normalized fibers with coverage (WLMAP), divided by the number of those
        fibers.

        The average spectrum is smoothed with a cubic spline with knots every
        `knot_step` Angstrom, fitted with weights proportional to the number of
        fibers with coverage, and only where at least `min_coverage` times the
        largest number of fibers have coverage; beyond, the spline keeps its
        value at the border. With `knot_step` <= 0, the average spectrum is
        smoothed with a spline of degree 5 with the default smoothing factor.

        A Savitzky-Golay filter (`window`, `degree`) of the average spectrum is
        only plotted in the intermediate results.
        """
        from scipy.signal import savgol_filter
        from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline

        # Bad fibers
        fp_conf = FocalPlaneConf.from_img(rss_wl)
        good_idxs_mask = numpy.ones((fp_conf.nfibers,), dtype="bool")
        good_idxs_mask[[fibid - 1 for fibid in fp_conf.invalid_fibers()]] = False

        data0 = rss_wl[0].data
        # This extension was created by WLcalibrator
        wlmap = rss_wl["WLMAP"].data > 0

        # Level of each fiber
        level = numpy.median(data0[:, col1:col2], axis=1)
        lit = good_idxs_mask & (level > 0)
        ref_level = numpy.median(level[lit]) if lit.any() else 0.0
        valid_mask = lit & (level > min_level * ref_level)
        faint = numpy.flatnonzero(lit & ~valid_mask) + 1
        if faint.size:
            self.logger.warning("fibers too faint for the average spectrum: %s", faint.tolist())

        valid_level = level[valid_mask]
        covered = wlmap[valid_mask]
        with numpy.errstate(invalid="ignore", divide="ignore"):
            data_good = numpy.where(covered, data0[valid_mask] / valid_level[:, numpy.newaxis], 0.0)
        data_good[~numpy.isfinite(data_good)] = 0.0

        # number of fibers in the average of each wavelength
        mm = covered.sum(axis=0)
        mask_noinfo = mm < 1
        collapse = numpy.where(mask_noinfo, 1.0, data_good.sum(axis=0) / numpy.maximum(mm, 1))
        # Smoothing works bad very near the border (overshooting)
        collapse_smooth = savgol_filter(collapse, window, degree)
        collapse_smooth[mask_noinfo] = 1.0

        xx = numpy.arange(collapse.shape[0])
        collapse_smooth_s = numpy.ones_like(collapse)
        if knot_step > 0:
            cdelt = abs(rss_wl[0].header.get("CDELT1", 1.0))
            coverage = mm / max(mm.max(), 1)
            fit = coverage >= min_coverage
            x_fit = xx[fit]
            nknots = max(int((x_fit[-1] - x_fit[0]) * cdelt / knot_step), 1)
            knots = numpy.linspace(x_fit[0], x_fit[-1], nknots + 1)[1:-1]
            self.logger.debug("average spectrum smoothed with %d interior knots", len(knots))
            spline = LSQUnivariateSpline(x_fit, collapse[fit], knots, w=coverage[fit], k=3, ext=3)
            spline_label = f"spline, knots every {knot_step:g} A"
        else:
            spline = UnivariateSpline(xx[~mask_noinfo], collapse[~mask_noinfo], k=5)
            spline_label = "spline5"
        collapse_smooth_s[~mask_noinfo] = spline(xx[~mask_noinfo])

        if self.intermediate_results:
            from matplotlib.figure import Figure

            numpy.savetxt("collapse.txt", collapse)
            numpy.savetxt("mask_noinfo.txt", mask_noinfo)
            fig = Figure()
            ax = fig.subplots()
            ax.plot(xx, collapse, ".", label="collapsed")
            ax.plot(xx, collapse_smooth, "-", label=f"savgol{degree}")
            ax.plot(xx, collapse_smooth_s, "--", label=spline_label)
            ax.legend()
            fig.savefig("collapsed_smooth.png")

        # Divide each fiber in rss_wl by spectrum
        gmean = valid_level.mean()
        data1 = rss_wl[0].data / collapse_smooth_s
        data1 /= gmean
        # Fill values with ones to avoid NaNs
        data2 = numpy.where(wlmap, data1, 1.0)

        self.logger.warning("Copy all extensions for the moment")
        rss_wl2 = fits.HDUList([hdu.copy() for hdu in rss_wl])
        rss_wl2[0].data = data2
        return rss_wl2

    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------
        rinput : RecipeInput

        Returns
        -------
        RecipeResult

        """

        img = self.process_flat2d(rinput)
        self.save_intermediate_img(img, "reduced_image.fits")
        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.master_apertures, self.datamodel, offset=rinput.extraction_offset)
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.master_wlcalib, self.datamodel)
        flipcor = FlipLR()

        img = splitter1(img)
        flat2d = splitter1.out  # Copy before extraction
        img = calibrator_aper(img)
        img = splitter2(img)
        rss_base = splitter2.out  # Copy before el calibration
        self.logger.debug("Flip RSS left-rigtht, before WL calibration")
        img = flipcor.run(img)
        # Calibrate in WL
        rss_wl = calibrator_wl(img)
        self.save_intermediate_img(rss_wl, "reduced_rss.fits")

        # Obtain flat field
        self.logger.info("Normalize flat field")
        rss_wl2 = self.obtain_fiber_flat(rss_wl, window=rinput.smoothing_window, knot_step=rinput.smoothing_knot_step)
        rss_wl2[0].header = self.set_product_headers(rss_wl2[0].header)
        result = self.create_result(master_fiberflat=rss_wl2, reduced_image=flat2d, reduced_rss=rss_base)
        return result

    def set_product_headers(self, hdr):
        """Set metadata in the FITS headers of the main product."""
        hdr = self.set_base_headers(hdr)
        hdr["NUMTYPE"] = ("MasterFiberFlat", "Product type")
        return hdr
