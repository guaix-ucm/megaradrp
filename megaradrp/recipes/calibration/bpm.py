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

import logging

from cosme import cosmetics, ccdmask, compute_mask
from numina.core import Product, DataFrameType
from numina.core.requirements import ObservationResultRequirement

import numpy as np

from megaradrp.core import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement

import matplotlib.pyplot as plt


_logger = logging.getLogger('numina.recipes.megara')


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
        for cont in range(1,1+len(rinput.obresult.images)//2):
            del(copia.images[fin - cont])
            del(copia2.images[0])

        reduced1 = self.bias_process_common(copia, rinput.master_bias)
        reduced2 = self.bias_process_common(copia2, rinput.master_bias)

        fits.writeto('reduced10.fits', reduced1[0].data, clobber=True)
        fits.writeto('reduced20.fits', reduced2[0].data, clobber=True)

        mask = np.zeros(reduced1[0].data.shape, dtype='int')
        bpm = cosmetics(reduced1[0].data, reduced2[1].data, mask)

        mask = np.zeros(reduced1[0].data.shape, dtype='int')
        bpm2 = ccdmask(reduced1[0].data, reduced2[0].data, mask, mode='full')

        bpm3 = compute_mask(reduced1[0].data, reduced2[0].data)

        fits.writeto('bpm3.fits', bpm3.astype(int), clobber=True)

        fits.writeto('maskresult1.fits', bpm[1], clobber=True)
        fits.writeto('result10.fits', bpm[0], clobber=True)

        fits.writeto('maskresult2.fits', bpm2[1], clobber=True)
        fits.writeto('result20.fits', bpm2[0], clobber=True)


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

    def run2(self, rinput):
        import matplotlib.pyplot as plt

        # Real detector
        SHAPE = (2056 * 2, 2048 * 2)
        PSCAN = 50
        OSCAN = 50

        shape_number = 500

        SHAPE = (shape_number, shape_number)
        # PSCAN = 4
        # OSCAN = 3

        # eq = 0.8 * np.ones(SHAPE)
        # eq[50,50:70] = 0.0
        # fits.writeto('eq.fits', eq, clobber=True)

        eq = np.random.normal(loc=0.80, scale=0.01, size=(shape_number,shape_number))
        eq = np.clip(eq, 0.0, 1.0)
        eq[50,50:70] = 0.0
        fits.writeto('eq.fits', eq, clobber=True)
        # eq = fits.getdata('eq.fits')

        # detector = Detector(eq=eq, gain=1.2, bias=100.0, ron=2.0)

        source = 5000
        number = 100

        accum = np.zeros(SHAPE)
        for frame in rinput.obresult.images:
            plt.imshow(frame, cmap='gray')
            plt.colorbar()

            # This is like bias subtraction
            # and overscan and trim
            cut = np.zeros_like(accum)
            cut[0:shape_number//2,0:shape_number] = frame[0:shape_number//2,PSCAN:shape_number+PSCAN]
            cut[shape_number//2:shape_number,0:shape_number] = frame[shape_number//2+2*OSCAN:shape_number+2*OSCAN,PSCAN:shape_number+PSCAN]
            plt.imshow(cut, cmap='gray')
            cut -= 100
            accum += cut

        lista_fin = accum / number

        mask = np.zeros(lista_fin.shape, dtype='int')
        result = ccdmask(lista_fin,None, mask)