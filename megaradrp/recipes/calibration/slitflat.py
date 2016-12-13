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

from __future__ import division, print_function

import uuid

import numpy
from scipy.ndimage.filters import median_filter
from astropy.io import fits
from numina.array import combine
from numina.core import Product, Parameter

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterSlitFlat
import megaradrp.requirements as reqs


class SlitFlatRecipe(MegaraBaseRecipe):
    """Process SLIT_FLAT images and create MASTER_SLIT_FLAT."""

    # Requirements
    master_bpm = reqs.MasterBiasRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    window_length_x = Parameter(301, 'Savitzky-Golay length of the filter window OX')
    window_length_y = Parameter(31, 'Savitzky-Golay length of the filter window OY')
    polyorder = Parameter(3, 'Savitzky-Golay order of the polynomial used to fit the samples')
    median_window_length = Parameter(31, 'Median window width')

    # Products
    master_slitflat = Product(MasterSlitFlat)

    def run(self, rinput):
        from scipy.signal import savgol_filter

        self.logger.info('starting slit flat reduction')

        flow = self.init_filters(rinput, rinput.obresult.configuration)
        reduced = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = reduced[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(reduced, 'reduced.fits')

        self.logger.debug("Compute median")
        if rinput.median_window_length:
            archivo_mediana = median_filter(reduced[0].data, (1, rinput.median_window_length))
        else:
            archivo_mediana = reduced[0].data

        self.logger.debug("Compute Savitzky-Golay X filter (%d, %d)", rinput.window_length_x, rinput.polyorder)
        result = savgol_filter(archivo_mediana, rinput.window_length_x, rinput.polyorder, axis=1)

        if rinput.window_length_y:
            self.logger.debug("Compute Savitzky-Golay Y filter (%d, %d)", rinput.window_length_y, rinput.polyorder)
            result = savgol_filter(result, rinput.window_length_y, rinput.polyorder, axis=0)

        qe = reduced[0].data / result

        self.logger.debug('Filtering INF/NAN in result')
        qe[numpy.isinf(qe)] = 1.0
        qe[numpy.isnan(qe)] = 1.0

        hdu = fits.PrimaryHDU(qe, header=reduced[0].header)
        hdu.header['UUID'] = uuid.uuid1().hex
        hdu.header['OBJECT'] = 'MASTER SLITFLAT'
        # hdu.header['IMAGETYP'] = 'SLIT_FLAT'

        master_slitflat = fits.HDUList([hdu])

        self.logger.info('end slit flat recipe')
        return self.create_result(master_slitflat=master_slitflat)
