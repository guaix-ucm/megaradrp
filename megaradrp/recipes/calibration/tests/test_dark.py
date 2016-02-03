from tempfile import mkdtemp
import numpy as np
import astropy.io.fits as fits
import pytest
import shutil

from megaradrp.tests.simulation import simulate_dark, simulate_dark_fits
from megaradrp.tests.simulation import ReadParams, MegaraDetectorSat, MegaraImageFactory
from numina.core import DataFrame, ObservationResult

from megaradrp.recipes.calibration.tests.test_bpm_common import generate_bias
from megaradrp.recipes.calibration.dark import DarkRecipe


def test_dark():
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    # ron = 2.0
    ron = 0.001
    gain = 1.0
    bias = 1000.0

    eq = 0.8 * np.ones(DSHAPE)

    temporary_path = mkdtemp()
    print ('Path: %s' %temporary_path)
    fits.writeto('%s/eq.fits' % temporary_path, eq, clobber=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    number = 10

    factory = MegaraImageFactory()
    fs = [simulate_dark_fits(factory, detector, exposure=3600) for i in range(number)]

    for aux in range(len(fs)):
        fs[aux].writeto('%s/dark_%s.fits' % (temporary_path, aux),clobber=True)

    master_bias = generate_bias(detector, number, temporary_path)
    master_bias_data = master_bias.master_bias.frame[0].data

    fits.writeto('%s/master_bias_data0.fits' % temporary_path,
                 master_bias_data, clobber=True)  # Master Bias

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    names = []

    for aux in range(number):
        names.append('%s/dark_%s.fits' % (temporary_path, aux))
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = DarkRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(temporary_path + '/master_bias_data0.fits').name))
    aux = recipe.run(ri)

    fits.writeto('%s/master_dark.fits' % temporary_path, aux.master_dark.frame[0].data, clobber=True)
    truncate_data = np.around(aux.master_dark.frame[0].data, decimals=2)
    shutil.rmtree(temporary_path)
    assert np.all(truncate_data == np.zeros(truncate_data.shape))

if __name__ == "__main__":
    test_dark()