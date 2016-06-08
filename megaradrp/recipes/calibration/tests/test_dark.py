import shutil
from tempfile import mkdtemp

import astropy.io.fits as fits
import numpy as np
from numina.core import DataFrame, ObservationResult

from megaradrp.recipes.calibration.dark import DarkRecipe
from megaradrp.recipes.calibration.tests.test_bpm_common import generate_bias
from megaradrp.simulation.factory import MegaraImageFactory
from megaradrp.simulation.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.actions import simulate_dark_fits
from megaradrp.core.insconf import MegaraInstrumentConfiguration

def test_dark():
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    # ron = 2.0
    ron = 0.001
    gain = 1.0
    bias = 1000.0

    qe = 0.8 * np.ones(DSHAPE)

    temporary_path = mkdtemp()
    fits.writeto('%s/eq.fits' % temporary_path, qe, clobber=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat('megara_test_detector', DSHAPE, OSCAN, PSCAN, qe=qe,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    number = 10

    factory = MegaraImageFactory()
    fs = simulate_dark_fits(factory, detector, exposure=3600, repeat=number)

    for idx, aux in enumerate(fs):
        aux.writeto('%s/dark_%s.fits' % (temporary_path, idx),clobber=True)

    master_bias = generate_bias(detector, number, temporary_path)
    master_bias_data = master_bias.master_bias.frame[0].data

    fits.writeto('%s/master_bias_data0.fits' % temporary_path,
                 master_bias_data, clobber=True)  # Master Bias

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.configuration = MegaraInstrumentConfiguration('configuration', {
        'trim1': [[0, 2056], [50, 4146]],
        'trim2': [[2156, 4212], [50, 4146]],
        'bng': [1, 1],
        'overscan1': [[0, 2056], [4149, 4196]],
        'overscan2': [[2156, 4212], [0, 50]],
        'prescan1': [[0, 2056], [0, 50]],
        'prescan2': [[2156, 4212], [4145, 4196]],
        'middle1': [[2056, 2106], [50, 4146]],
        'middle2': [[2106, 2156], [50, 4146]]})
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