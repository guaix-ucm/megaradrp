#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Tests for the bpm mode recipe module."""
import shutil
from tempfile import mkdtemp

import astropy.io.fits as fits
import numpy as np

from numina.core import DataFrame, ObservationResult
import numina.instrument.assembly as asb

from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe
from megaradrp.recipes.calibration.tests.test_bpm_common import generate_bias
from megaradrp.instrument.components.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.actions import simulate_flat


# @pytest.mark.remote_data
def test_bpm():
    number = 5
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    ron = 2.0
    gain = 1.0
    bias = 1000.0

    qe = 0.8 * np.ones(DSHAPE)
    qe[5:6, 0:170] = 0.0
    config_uuid = '4fd05b24-2ed9-457b-b563-a3c618bb1d4c'
    date_obs = '2017-11-09T11:00:00.0'
    temporary_path = mkdtemp()

    fits.writeto('{}/eq.fits'.format(temporary_path), qe, overwrite=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat('megara_test_detector', DSHAPE, OSCAN, PSCAN,
                                 qe=qe,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    source2 = 1.0

    fs = [simulate_flat(detector, exposure=1.0, source=5000 * source2) for i in
          range(number)]
    fs2 = [simulate_flat(detector, exposure=1.0, source=40000 * source2) for i
           in range(number)]

    header = fits.Header()
    header['DATE-OBS'] = date_obs
    header['INSCONF'] = config_uuid
    header['INSTRUME'] = 'MEGARA'
    header['VPH'] = 'LR-U'
    header['INSMODE'] = 'MOS'
    for aux in range(len(fs)):
        fits.writeto('{}/flat_{}.fits'.format(temporary_path, aux), fs[aux],
                     header=header,
                     overwrite=True)
        fits.writeto('{}/flat_{}.fits'.format(temporary_path, aux + number), fs2[aux],
                     header=header,
                     overwrite=True)

    result = generate_bias(detector, number, temporary_path)
    result.master_bias.frame.writeto(
        '{}/master_bias_data0.fits'.format(temporary_path),
        overwrite=True
    )

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'MegaraBiasImage'
    pkg_paths = ['megaradrp.instrument.configs']
    store = asb.load_paths_store(pkg_paths)
    insmodel = asb.assembly_instrument(store, config_uuid, date_obs, by_key='uuid')
    insmodel.configure_with_header(header)
    ob.configuration = insmodel

    names = []
    for aux in range(number * 2):
        names.append('{}/flat_{}.fits'.format(temporary_path, aux))
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(temporary_path + '/master_bias_data0.fits').name))
    aux = recipe.run(ri)
    aux.master_bpm.frame.writeto('{}/master_bpm.fits'.format(temporary_path), overwrite=True)
    shutil.rmtree(temporary_path)


if __name__ == "__main__":
    test_bpm()
