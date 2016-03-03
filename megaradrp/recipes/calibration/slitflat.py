#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

import logging
import numpy as np

from astropy.io import fits

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.products import MasterSlitFlat
from megaradrp.requirements import MasterFiberFlatRequirement
from numina.core import Product, Parameter
from scipy.ndimage.filters import median_filter
from scipy.signal import savgol_filter

_logger = logging.getLogger('numina.recipes.megara')

class SlitFlatRecipe(MegaraBaseRecipe):
    """Process SLIT_FLAT images and create MASTER_SLIT_FLAT."""

    # Requirements
    master_fiberflat = MasterFiberFlatRequirement()
    window_length_x = Parameter(301, 'Savitzky-Golay length of the filter window OX')
    window_length_y = Parameter(31, 'Savitzky-Golay length of the filter window OY')
    polyorder = Parameter(3, 'Savitzky-Golay order of the polynomial used to fit the samples')
    median_window_length = Parameter(31, 'Median window width')

    # Products
    master_slitflat = Product(MasterSlitFlat)

    def __init__(self):
        super(SlitFlatRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        _logger.info('Slit Flat')

        with rinput.master_fiberflat.open() as hdul:
            master_fiberflat = hdul[0].data.copy()

        master_fiberflat[0:2055,:] = master_fiberflat[0:2055,:]-1000.0
        master_fiberflat[2056:4111,:] = master_fiberflat[2056:4111,:]-1005.0

        if rinput.median_window_length:
            archivo_mediana = median_filter(master_fiberflat, (1,rinput.median_window_length))
        else:
            archivo_mediana = master_fiberflat

        result = savgol_filter(archivo_mediana, rinput.window_length_x, rinput.polyorder, axis=1)

        if rinput.window_length_y:
            result = savgol_filter(result, rinput.window_length_y, rinput.polyorder, axis=0)

        qe = master_fiberflat/result

        qe[np.isinf(qe)] = 1.0
        qe[np.isnan(qe)] = 1.0

        hdu = fits.PrimaryHDU(qe)
        reduced = fits.HDUList([hdu])

        return self.create_result(master_slitflat=reduced)
