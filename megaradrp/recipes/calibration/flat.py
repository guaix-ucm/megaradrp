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
from megaradrp.processing.datamodel import MegaraDataModel
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator


from numina.flow.processing import Corrector
class Splitter(Corrector):
    """Make a copy of its input"""
    def __init__(self):
        super(Splitter, self).__init__()
        self.out = None

    def run(self, img):
        self.out = self.copy_img(img)
        return img

    def copy_img(self, img):
        return fits.HDUList([hdu.copy() for hdu in img])


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

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(version="0.1.0")

    def process_flat2d(self, rinput):
        flow = self.init_filters(rinput, rinput.obresult.configuration.values)
        final_image = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = final_image[0].header
        self.set_base_headers(hdr)
        return final_image

    def obtain_fiber_flat(self, rss_wl, wlcalib):
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
        import matplotlib.pyplot as plt
        # FIXME: skip incorrect traces
        data_good = data0[valid_mask] / col_good_mean[:, numpy.newaxis]
        data_good[numpy.isnan(data_good)] = 0.0
        plt.imshow(data_good)
        plt.show()
        collapse = numpy.mean(data_good, axis=0)
        import matplotlib.pyplot as plt
        plt.plot(collapse)
        plt.show()
        # Filter collapse to smooth it
        # Divide each fiber in rss_wl by the smoothed spectrum

        gmean = col_good_mean.mean()
        data1 = rss_wl[0].data / collapse
        data1 /= gmean

        # import matplotlib.pyplot as plt
        # for m in data_good:
        #     plt.plot(m, 'k.')
        # plt.plot(collapse, 'b-')
        # plt.show()
        # plt.imshow(data1)
        # plt.show()

        rss_wl2 = fits.HDUList([hdu.copy() for hdu in rss_wl])
        rss_wl2[0].data = data1
        return rss_wl2
        # self.logger.info('Compute mean and resampling again')

    def run(self, rinput):

        img = self.process_flat2d(rinput)


        datamodel = MegaraDataModel()
        splitter1 = Splitter()
        calibrator_aper = ApertureExtractor(rinput.tracemap, datamodel)
        splitter2 = Splitter()
        calibrator_wl = WavelengthCalibrator(rinput.wlcalib, datamodel)

        img = splitter1(img)
        flat2d = splitter1.out # Copy before extraction
        img = calibrator_aper(img)
        img = splitter2(img)
        rss_base = splitter2.out # A before el calibration

        # FIXME: Flip L-R image before calibrating WL
        # Eventually this should not be necessary
        import numpy
        self.logger.debug('Flip RSS left-rigtht, before WL calibration')
        img[0].data = numpy.fliplr(img[0].data)

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

