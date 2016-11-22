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

"""Fiber flat calibration Recipe for Megara"""

from __future__ import division, print_function

import numpy
from astropy.io import fits

from numina.core import Product, Requirement
from numina.flow import SerialFlow

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterFiberFlat
from megaradrp.products import WavelengthCalibration
import megaradrp.requirements as reqs
from megaradrp.types import ProcessedRSS, ProcessedFrame

# Flat 2D
from megaradrp.processing.combine import basic_processing_with_combination
from numina.array import combine
# Create RSS
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR

class FiberFlatRecipe(MegaraBaseRecipe):
    """Process FIBER_FLAT images and create MASTER_FIBER_FLAT."""

    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    wlcalib = Requirement(WavelengthCalibration, 'Wavelength calibration table')
    # master_weights = Requirement(MasterWeights, 'Set of files')
    tracemap = reqs.MasterTraceMapRequirement()

    # Products
    fiberflat_frame = Product(ProcessedFrame)
    fiberflat_rss = Product(ProcessedRSS)
    master_fiberflat = Product(MasterFiberFlat)

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        final_image = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def obtain_fiber_flat(self, rss_wl, wlcalib):
        from scipy.signal import savgol_filter

        # Bad fibers, join:
        bad_fibers = wlcalib.missing_fibers
        bad_fibers.extend(wlcalib.error_fitting)
        # print(bad_fibers)
        bad_idxs = [fibid - 1 for fibid in bad_fibers]
        # print(bad_idxs)

        good_idxs_mask = numpy.ones((wlcalib.total_fibers,), dtype='bool')
        good_idxs_mask[bad_idxs] = False

        # Collapse all fiber spectrum
        xcol = slice(1900, 2100)  # FIXME: hardcoded

        data0 = rss_wl[0].data

        col_mean = data0[:, xcol].mean(axis=1)
        # Filter positive values and valid fibers
        col_mean_pos = (col_mean > 0)
        valid_mask = col_mean_pos & good_idxs_mask

        col_good_mean = col_mean[valid_mask]

        # FIXME: skip incorrect traces
        data_good = data0[valid_mask] / col_good_mean[:, numpy.newaxis]
        data_good[numpy.isnan(data_good)] = 0.0

        # Crappy way
        # This extension was created by WLcalibrator
        wlmap = rss_wl['WLMAP'].data
        mm = numpy.sum(wlmap, axis=0)

        # Filter collapse to smooth it
        collapse = numpy.sum(data_good, axis=0) / mm
        # FIXME: Savitsky-Golay filter, window 31, pol degree 3
        collapse_smooth = savgol_filter(collapse, 31, 3)

        # Divide each fiber in rss_wl by  spectrum
        gmean = col_good_mean.mean()
        data1 = rss_wl[0].data / collapse_smooth
        data1 /= gmean
        # Fill values with ones to avoid NaNs
        data2 = numpy.where(wlmap > 0, data1, 1.0)

        self.logger.warning("Copy all extensions for the moment")
        rss_wl2 = fits.HDUList([hdu.copy() for hdu in rss_wl])
        rss_wl2[0].data = data2
        return rss_wl2

    def run(self, rinput):

        img = self.process_flat2d(rinput)
        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.tracemap, self.datamodel)
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.wlcalib, self.datamodel)
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

        # Obtain flat field
        self.logger.info('Normalize flat field')
        rss_wl2 = self.obtain_fiber_flat(rss_wl, rinput.wlcalib)
        result = self.create_result(
            master_fiberflat=rss_wl2,
            fiberflat_frame=flat2d,
            fiberflat_rss=rss_base
        )
        return result
