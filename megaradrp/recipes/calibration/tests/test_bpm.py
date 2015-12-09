#
# Copyright 2015 Universidad Complutense de Madrid
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

"""Tests for the arc mode recipe module."""
import astropy.io.fits as fits
import numpy as np
import pytest

from numina.array.cosmetics import cosmetics, ccdmask
from numina.core import DataFrame, ObservationResult
from numina.tests.testcache import download_cache

from megaradrp.tests.simulation import MegaraDetector as Detector
from megaradrp.tests.simulation import simulate_flat, simulate_bias
from megaradrp.tests.simulation import ReadParams, MegaraImageFactory, MegaraDetectorSat
from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe
from megaradrp.recipes.calibration.bias import BiasRecipe

from megaradrp.requirements import MasterBiasRequirement

import matplotlib.pyplot as plt
from sklearn import svm



def __test_generate_flats():
    PSCAN = 50

    # DSHAPE = (2056 * 2, 2048 * 2)  # Original
    DSHAPE = (500 * 2, 500 * 2)
    OSCAN = 0

    ron = 2.0
    gain = 1.0
    bias = 1000.0

    # eq = np.random.normal(loc=0.80, scale=0.01, size=SHAPE)
    # eq = np.clip(eq, 0.0, 1.0)
    eq = 0.8 * np.ones(DSHAPE)
    eq[50:55,250:770] = 0.0
    fits.writeto('eq.fits', eq, clobber=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq, dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2, bins='11')

    # Non uniform illumination
    x = np.linspace(-1.0, 1.0, num=eq.shape[1])
    y = np.linspace(-1.0, 1.0, num=eq.shape[0])
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = 1-0.3*(xx**2 + yy**2)
    source1 = z

    number = 30
    fs = [simulate_flat(detector, exposure=1.0, source=5000*source1).astype('int')-1000 for i in range(number)]
    # fs2 = [simulate_flat(detector, exposure=1.0, source=15000*source1).astype('int')-1000 for i in range(number)]
    aux = np.array(fs)
    # aux2 = np.array(fs2)
    med_arr = np.median(aux, axis=0)
    # med_arr2 = np.median(aux2, axis=0)
    var_arr = np.var(aux, axis=0)
    # var_arr2 = np.var(aux2, axis=0)
    z = np.polyfit(med_arr.ravel(), var_arr.ravel(), 1)

    import matplotlib.pyplot as plt

    print z[0], z[1]

    # total = np.arange(12500)
    # plt.plot(med_arr.ravel(), var_arr.ravel(),'.',alpha=0.01)
    # plt.plot(med_arr2.ravel(), var_arr2.ravel(),'.',alpha=0.01)
    # plt.plot(total, z[0]*total+z[1])
    # plt.show()

    bad_pixels = (np.array(med_arr.ravel())<=1.0).astype('int')







def generate_bias(detector, number):


    fs = [simulate_bias(detector) for i in range(number)]
    # eq[50:55,250:770] = 0.0
    # mask2 = np.ones_like(fs[0])
    # mask2[50:55,250:770] = 0
    # fs = [f*mask2 for f in fs]
    # #
    for aux in range(len(fs)):
        fits.writeto('bias_%s.fits'%aux, fs[aux], clobber=True)

    fs = ["bias_%s.fits"%i for i in range(number)]

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=f) for f in fs]

    recipe = BiasRecipe()
    ri = recipe.create_input(obresult=ob)
    return recipe.run(ri)

def test_bpm():
    # PSCAN = 50
    # OSCAN = 50
    # shape_number = 1000
    number = 5
    # ron = 2.0
    # gain = 1.0
    # bias = 1000.0
    #
    # OSCAN = 100
    PSCAN = 50
    # H_X_DIM = 2048 * 2
    # H_Y_DIM = 2056 * 2
    #
    # SHAPE1 = (H_X_DIM, H_Y_DIM)
    # SHAPE = (H_Y_DIM,H_X_DIM)

    DSHAPE = (2056 * 2, 2048 * 2)
    # PSCAN = 0
    OSCAN = 50

    BINR = 1
    BINC = 1

    SHAPE = DSHAPE[0] // BINR, DSHAPE[1] // BINC
    HSHAPE = SHAPE[0] // 2, DSHAPE[1] // 2

    ron = 2.0
    gain = 1.0
    bias = 1000.0

    # eq = np.random.normal(loc=0.80, scale=0.01, size=SHAPE)
    # eq = np.clip(eq, 0.0, 1.0)
    eq = 0.8 * np.ones(DSHAPE)
    eq[50:55,250:770] = 0.0
    fits.writeto('eq.fits', eq, clobber=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    # factory = MegaraImageFactory()

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq, dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2, bins='11')

    source2 = 1.0

    # raw12 = np.zeros_like(eq)
    # raw22 = np.zeros_like(eq)

    fs = [simulate_flat(detector, exposure=1.0, source=5000*source2) for i in range(number)]
    fs2 = [simulate_flat(detector, exposure=1.0, source=40000*source2) for i in range(number)]

    for aux in range(len(fs)):
        fits.writeto('flat_%s.fits'%aux, fs[aux], clobber=True)
        fits.writeto('flat_%s.fits'%(aux+number), fs2[aux], clobber=True)

    master_bias = generate_bias(detector, number)
    master_bias_data = master_bias.biasframe.frame[0].data
    # plt.imshow(master_bias.biasframe.frame[2].data, cmap='gray')
    fits.writeto('master_bias_data0.fits', master_bias_data, clobber=True)  # Master Bias
    fits.writeto('master_bias_data1.fits', master_bias.biasframe.frame[1].data, clobber=True)  # Master Bias
    fits.writeto('master_bias_data2.fits', master_bias.biasframe.frame[2].data, clobber=True)  # Master Bias

    # source = 5000

    # detector = Detector(shape=SHAPE, oscan=OSCAN, pscan=PSCAN, eq=eq, gain=1.2, bias=100.0, ron=2.0)
    # fs = [simulate_flat(detector, exposure=10.0, source=source) for i in range(number)]
    # mask2 = np.ones_like(fs[0])
    # mask2[0:500,300] = 0
    # fs = [f*mask2 for f in fs]
    # cont = 0
    # for aux in range(len(fs)):
    #     fits.writeto('flat_%s.fits'%aux, fs[aux], clobber=True)
    #     cont +=1
    # del(fs)

    # fs = [simulate_flat(detector, exposure=15.0, source=source) for i in range(number)]
    # mask2 = np.ones_like(fs[0])
    # mask2[0:500,300] = 0
    # fs = [f*mask2 for f in fs]
    #
    # for aux in range(len(fs)):
    #     fits.writeto('flat_%s.fits'%cont, fs[aux], clobber=True)
    #     cont +=1

    # n_imagenes = 5
    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    names = []
    for aux in range(number*2):
        names.append('flat_%s.fits' %(aux))
    # names = ['flat_%s.fits' %cont for cont in range(number*2)]
    ob.frames = [DataFrame(filename=file(nombre).name) for nombre in names]

    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(filename=file('master_bias_data0.fits').name))
    recipe.run(ri)




def _test_bpm(drpmocker):

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

    detector = Detector(eq=eq, gain=1.2, bias=100.0, ron=2.0)

    source = 5000
    number = 100

    accum = np.zeros(SHAPE)
    lista_fin = []
    for exposure_time in [10.0,15.0]:
        for i in range(number):

            raw = simulate_flat(detector, exposure=exposure_time, source=source)
            fits.writeto('raw.fits', raw, clobber=True)

            print('imag', i, raw.mean())
            # This is like bias subtraction
            # and overscan and trim
            cut = np.zeros_like(accum)
            cut[0:shape_number//2,0:shape_number] = raw[0:shape_number//2,PSCAN:shape_number+PSCAN]
            cut[shape_number//2:shape_number,0:shape_number] = raw[shape_number//2+2*OSCAN:shape_number+2*OSCAN,PSCAN:shape_number+PSCAN]

            cut -= 100
            accum += cut
            fits.writeto('cut.fits', cut, clobber=True)

        lista_fin.append(accum / number)
        accum = np.zeros(SHAPE)



    # dark_median_pica = np.median(np.array(lista2),axis=0)
    # fits.writeto('dark_median_pica.fits', dark_median_pica, clobber=True)


    fits.writeto('lista_fin0.fits', lista_fin[0], clobber=True)
    fits.writeto('lista_fin1.fits', lista_fin[1], clobber=True)
    fits.writeto('ratio.fits', lista_fin[1] / lista_fin[0], clobber=True)



    mask = np.zeros(lista_fin[0].shape, dtype='int')
    fits.writeto('mask_original.fits', mask, clobber=True)

    result = cosmetics(lista_fin[0], lista_fin[1], mask)

    fits.writeto('maskresult1.fits', mask, clobber=True)
    fits.writeto('result10.fits', result[0], clobber=True)
    fits.writeto('result11.fits', result[1], clobber=True)

    mask = np.zeros(lista_fin[0].shape, dtype='int')
    result2 = ccdmask(lista_fin[0],None, mask)
    fits.writeto('maskresult2.fits', mask, clobber=True)
    fits.writeto('result20.fits', result2[0], clobber=True)
    fits.writeto('result21.fits', result2[1], clobber=True)



if __name__ == "__main__":
    # test_bpm()
    test_generate_flats()