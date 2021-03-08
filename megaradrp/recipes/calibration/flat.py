#
# Copyright 2011-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Fiber flat calibration Recipe for Megara"""

from __future__ import division, print_function

import numpy
from astropy.io import fits
import matplotlib.pyplot as plt
from numina.core import Result, Parameter
import numina.exceptions

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.ntypes import MasterFiberFlat
import megaradrp.requirements as reqs
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame

# Flat 2D
from megaradrp.processing.combine import basic_processing_with_combination
from numina.array import combine
# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR


def _smoothing_window_check(val):
    """Check 'smoothing_window' """
    #raise ValueError("something")
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
    method = Parameter(
        'median',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_apertures = reqs.MasterAperturesRequirement(alias='master_traces')
    smoothing_window = Parameter(31, 'Window for smoothing (must be odd)',
                                 validator=_smoothing_window_check
                                 )
    extraction_offset = Parameter([0.0], 'Offset traces for extraction', accept_scalar=True)
    master_wlcalib = reqs.WavelengthCalibrationRequirement()

    # Results
    reduced_image = Result(ProcessedFrame)
    reduced_rss = Result(ProcessedRSS)
    master_fiberflat = Result(MasterFiberFlat)

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        fmethod = getattr(combine, rinput.method)
        final_image = basic_processing_with_combination(
            rinput, flow, method=fmethod, method_kwargs=rinput.method_kwargs,
        )
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def obtain_fiber_flat(self, rss_wl, col1=1900, col2=2100, window=31, degree=3):
        from scipy.signal import savgol_filter
        from scipy.interpolate import UnivariateSpline

        # Bad fibers
        fp_conf = FocalPlaneConf.from_img(rss_wl)
        bad_fibers = fp_conf.invalid_fibers()
        bad_idxs = [fibid - 1 for fibid in bad_fibers]
        # print(bad_idxs)

        good_idxs_mask = numpy.ones((fp_conf.nfibers,), dtype='bool')
        good_idxs_mask[bad_idxs] = False

        # Collapse all fiber spectrum
        xcol = slice(col1, col2)

        data0 = rss_wl[0].data

        col_mean = data0[:, xcol].mean(axis=1)
        # Filter positive values and valid fibers
        col_mean_pos = (col_mean > 0)
        valid_mask = col_mean_pos & good_idxs_mask

        col_good_mean = col_mean[valid_mask]

        data_good = data0[valid_mask] / col_good_mean[:, numpy.newaxis]
        data_good[numpy.isnan(data_good)] = 0.0

        # This extension was created by WLcalibrator
        wlmap = rss_wl['WLMAP'].data
        mm = numpy.sum(wlmap, axis=0)
        # The information is also in the keywords
        # FIBxxxS1, FIBxxxS2
        # skip 0 in divisions
        mask_noinfo = mm < 1
        mm[mask_noinfo] = 1
        # Filter collapse to smooth it
        collapse = numpy.sum(data_good, axis=0) / mm
        # Smoothing works bad very near the border (overshooting)
        collapse_smooth = savgol_filter(collapse, window, degree)
        collapse_smooth[mask_noinfo] = 1.0

        xx = numpy.arange(collapse.shape[0])
        xx_v = xx[~mask_noinfo]
        collapse_v = collapse[~mask_noinfo]

        # Spline degree (3 or 5)
        degree_s = 5
        interpol_univ = UnivariateSpline(xx_v, collapse_v, k=degree_s)
        collapse_smooth_s = collapse_smooth.copy()
        collapse_smooth_s[~mask_noinfo] = interpol_univ(xx_v)
        collapse_smooth_s[mask_noinfo] = 1.0

        if self.intermediate_results:
            numpy.savetxt('collapse.txt', collapse)
            numpy.savetxt('mask_noinfo.txt', mask_noinfo)
            fig, ax = plt.subplots()
            ax.plot(xx, collapse, '.', label='collapsed')
            ax.plot(xx, collapse_smooth, '-', label=f'savgol{degree}')
            ax.plot(xx, collapse_smooth_s, '--', label=f'spline{degree_s}')
            ax.legend()
            plt.savefig('collapsed_smooth.png')
            plt.close()

        # Divide each fiber in rss_wl by spectrum
        gmean = col_good_mean.mean()
        data1 = rss_wl[0].data / collapse_smooth_s
        data1 /= gmean
        # Fill values with ones to avoid NaNs
        data2 = numpy.where(wlmap > 0, data1, 1.0)

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
        self.save_intermediate_img(img, 'reduced_image.fits')
        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(
            rinput.master_apertures,
            self.datamodel,
            offset=rinput.extraction_offset
        )
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.master_wlcalib, self.datamodel)
        flipcor = FlipLR()

        img = splitter1(img)
        flat2d = splitter1.out # Copy before extraction
        img = calibrator_aper(img)
        img = splitter2(img)
        rss_base = splitter2.out # Copy before el calibration
        self.logger.debug('Flip RSS left-rigtht, before WL calibration')
        img = flipcor.run(img)
        # Calibrate in WL
        rss_wl = calibrator_wl(img)
        self.save_intermediate_img(rss_wl, 'reduced_rss.fits')

        # Obtain flat field
        self.logger.info('Normalize flat field')
        rss_wl2 = self.obtain_fiber_flat(rss_wl, window=rinput.smoothing_window)
        rss_wl2[0].header = self.set_base_headers(rss_wl2[0].header)
        result = self.create_result(
            master_fiberflat=rss_wl2,
            reduced_image=flat2d,
            reduced_rss=rss_base
        )
        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(FiberFlatRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterFiberFlat', 'Product type')
        return hdr