#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import shutil
from tempfile import mkdtemp

import astropy.io.fits as fits
import numpy
from numina.core import DataFrame, ObservationResult
import numina.instrument.assembly as asb

from megaradrp.recipes.calibration.dark import DarkRecipe
from megaradrp.simulation.factory import MegaraImageFactory
from megaradrp.instrument.components.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.actions import simulate_dark_fits


def test_dark():

    numpy.random.seed(422992983)

    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    # ron = 2.0
    ron = 0.001
    gain = 1.0
    bias = 1000.0
    dark = 3.0 # In 1 hour
    exptime = 3600.0
    dark_s = dark / exptime
    qe = 0.8 * numpy.ones(DSHAPE)
    config_uuid = '4fd05b24-2ed9-457b-b563-a3c618bb1d4c'
    date_obs = '2017-11-09T11:00:00.0'
    temporary_path = mkdtemp()
    fits.writeto('{}/eq.fits'.format(temporary_path), qe, overwrite=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat('megara_test_detector', DSHAPE, OSCAN, PSCAN, qe=qe,
                                 dark=dark_s,
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    number = 3
    factory = MegaraImageFactory()
    fs = simulate_dark_fits(factory, detector, exposure=3600, repeat=number)

    for idx, aux in enumerate(fs):
        aux.writeto('{}/dark_{}.fits'.format(temporary_path, idx), overwrite=True)

    header = fits.Header()
    header['DATE-OBS'] = date_obs
    header['INSCONF'] = config_uuid
    header['INSTRUME'] = 'MEGARA'
    header['VPH'] = 'LR-U'
    header['INSMODE'] = 'MOS'
    master_bias_data = numpy.zeros(DSHAPE)
    master_bias_hdul = fits.HDUList(fits.PrimaryHDU(
        master_bias_data, header=header)
    )
    #master_bias_data = master_bias.master_bias.frame[0].data

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'MegaraDarkImage'

    pkg_paths = ['megaradrp.instrument.configs']
    store = asb.load_paths_store(pkg_paths)
    insmodel = asb.assembly_instrument(store, config_uuid, date_obs, by_key='uuid')
    insmodel.configure_with_header(header)
    ob.configuration = insmodel

    names = []
    for aux in range(number):
        names.append('{}/dark_{}.fits'.format(temporary_path, aux))
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    ob.configuration = insmodel
    recipe = DarkRecipe()
    ri = recipe.create_input(
        obresult=ob,
        master_bias=DataFrame(frame=master_bias_hdul),
    )
    aux = recipe.run(ri)

    mean_dark_value = aux.master_dark.frame[0].data.mean()

    shutil.rmtree(temporary_path)

    assert numpy.allclose(mean_dark_value, dark, atol=0, rtol=1e-1)


if __name__ == "__main__":
    test_dark()
