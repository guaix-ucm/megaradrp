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
#

"""Bad PIxel Mask (BPM) recipe"""

import astropy.io.fits as fits
import numpy as np
from numina.array.cosmetics import ccdmask
from numina.core import Product, DataFrameType
from numina.core.requirements import ObservationResultRequirement

from megaradrp.core import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement


class BadPixelsMaskRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()

    bpm_image = Product(DataFrameType)

    def __init__(self):
        super(BadPixelsMaskRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        import copy

        N = len(rinput.obresult.frames)
        obresult1 = copy.copy(rinput.obresult)
        obresult1.frames = rinput.obresult.frames[:N//2]
        obresult2 = copy.copy(rinput.obresult)
        obresult2.frames = rinput.obresult.frames[N//2:]

        reduced1 = self.bias_process_common(obresult1, rinput.master_bias)
        reduced2 = self.bias_process_common(obresult2, rinput.master_bias)

        mask = np.zeros(reduced1[0].data.shape, dtype='int')

        bpm = ccdmask(reduced1[0].data, reduced2[0].data, mask, mode='full')
        hdu = fits.PrimaryHDU(bpm)

        reduced = fits.HDUList([hdu])

        return self.create_result(bpm_image=reduced)
