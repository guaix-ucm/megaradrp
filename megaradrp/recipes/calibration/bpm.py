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
import os

from megaradrp.core import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement

from numina.array.cosmetics import cosmetics, ccdmask
from numina.core import Product, DataFrameType
from numina.core.requirements import ObservationResultRequirement

class BadPixelsMaskRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()

    bpm_image = Product(DataFrameType)

    def __init__(self):
        super(BadPixelsMaskRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):
        import copy
        copia = copy.deepcopy(rinput.obresult)
        copia2 = copy.deepcopy(rinput.obresult)
        fin = len(rinput.obresult.images)
        for cont in range(1, 1 + len(rinput.obresult.images) // 2):
            del (copia.images[fin - cont])
            del (copia2.images[0])

        temporary_path = os.path.dirname(
            os.path.abspath(copia.images[0].filename))

        reduced1 = self.bias_process_common(copia, rinput.master_bias)
        reduced2 = self.bias_process_common(copia2, rinput.master_bias)

        fits.writeto('%s/reduced10.fits' % temporary_path, reduced1[0].data,
                     clobber=True)
        fits.writeto('%s/reduced20.fits' % temporary_path, reduced2[0].data,
                     clobber=True)

        mask = np.zeros(reduced1[0].data.shape, dtype='int')
        bpm = cosmetics(reduced1[0].data, reduced2[0].data, mask)
        print '# bad pixels cosmetics', bpm.sum()

        mask = np.zeros(reduced1[0].data.shape, dtype='int')
        bpm2 = ccdmask(reduced1[0].data, reduced2[0].data, mask, mode='full')
        print '# bad pixels ccdmask', bpm2[1].sum()

        fits.writeto('%s/maskresult1.fits' % temporary_path, bpm.astype(int),
                     clobber=True)

        fits.writeto('%s/maskresult2.fits' % temporary_path, bpm2[1],
                     clobber=True)
        fits.writeto('%s/result20.fits' % temporary_path, bpm2[0],
                     clobber=True)


        # hdu, data = self.hdu_creation(rinput.obresult, {'bpm':bpm})
        #
        # hdr = hdu.header
        # # FIXME: this is incorrect in general
        # hdr = self.set_base_headers(hdr)
        # hdr['CCDMEAN'] = data[0].mean()
        #
        # varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        # num = fits.ImageHDU(data[2], name='MAP')
        #
        # reduced = fits.HDUList([hdu, varhdu, num])
        #
        # return self.create_result(bpm_image=reduced)
