#
# Copyright 2011-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division

import uuid

import numpy
from scipy.ndimage import median_filter
from astropy.io import fits
from numina.array import combine
from numina.core import Result

from megaradrp.ntypes import ProcessedFrame
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import MasterSlitFlat
import megaradrp.requirements as reqs


class SlitFlatRecipe(MegaraBaseRecipe):
    """Process SLIT_FLAT images and create MasterSlitFlat."""

    # Requirements
    master_bpm = reqs.MasterBPMRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    # window_length_x = Parameter(301, 'Savitzky-Golay length of the filter window OX')
    # window_length_y = Parameter(31, 'Savitzky-Golay length of the filter window OY')
    # polyorder = Parameter(3, 'Savitzky-Golay order of the polynomial used to fit the samples')
    # median_window_length = Parameter(31, 'Median window width')

    # Results
    reduced_image = Result(ProcessedFrame)
    master_slitflat = Result(MasterSlitFlat)

    def run(self, rinput):
        self.logger.info('starting slit flat reduction')

        flow = self.init_filters(rinput, rinput.obresult.configuration)
        reduced = basic_processing_with_combination(
            rinput, flow, method=combine.median)
        hdr = reduced[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(reduced, 'reduced_image.fits')

        # Using median filtering... In each channel
        l0 = reduced[0].data.shape[0]
        l1 = l0 // 2
        channel1 = (slice(None, l1), slice(None, None, None))
        channel2 = (slice(l1, None, None), slice(None, None, None))

        m_window1 = (11, 11)
        self.logger.debug('median filtering by channel %s', m_window1)
        median1 = numpy.zeros_like(reduced[0].data)
        median1[channel1] = median_filter(reduced[0].data[channel1], m_window1)
        median1[channel2] = median_filter(reduced[0].data[channel2], m_window1)

        self.save_intermediate_array(median1, 'median_image1.fits')

        qe1 = reduced[0].data / median1

        self.save_intermediate_array(qe1, 'qe_filter1.fits')

        m_window2 = (3, 21)
        self.logger.debug('median filtering %s', m_window2)

        median2 = median_filter(qe1, m_window2)
        self.save_intermediate_array(median2, 'median_image2.fits')

        qe2 = qe1 / median2
        self.save_intermediate_array(qe2, 'qe_filter2.fits')

        self.logger.debug('filtering Inf/NaN in result')
        qe2[numpy.isinf(qe2)] = 1.0
        qe2[numpy.isnan(qe2)] = 1.0

        hdu = fits.PrimaryHDU(qe2.astype('float32'), header=reduced[0].header)
        hdu.header['UUID'] = str(uuid.uuid1())

        master_slitflat = fits.HDUList([hdu])
        self.set_base_headers(master_slitflat[0].header)

        self.logger.info('end slit flat recipe')
        return self.create_result(master_slitflat=master_slitflat,
                                  reduced_image=reduced)

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(SlitFlatRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterSlitFlat', 'Product type')
        return hdr
