#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import os

import pytest

from numina.tests.testcache import download_cache
from numina.core import ObservationResult
from numina.core import DataFrame
from numina.instrument.assembly import assembly_instrument

from megaradrp.recipes.calibration.bias import BiasRecipe
from megaradrp.loader import load_drp


@pytest.mark.remote_data
def test_bias():

    BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/BIAS/%s'
    images = ['e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e']

    fs = [download_cache(BASE_URL % i) for i in images]

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'MegaraBiasImage'
    ob.frames = [DataFrame(filename=f.name) for f in fs]

    insdrp = load_drp()
    pipeline = insdrp.pipelines.get('default')
    recipe = pipeline.get_recipe_object(ob.mode)

    assert isinstance(recipe, BiasRecipe)
    # TODO: these should be created by the DAL

    import numina.instrument.assembly as asbl
    import numina.drps
    sys_drps = numina.drps.get_system_drps()
    com_store = asbl.load_panoply_store(sys_drps)

    key, date_obs, keyname = insdrp.select_profile(ob)
    ob.configuration = assembly_instrument(com_store, key, date_obs, by_key=keyname)

    ri = recipe.create_input(obresult=ob)

    result = recipe.run(ri)
    # assert result.qc >= QC.UNKNOWN

    # Checks on the image
    hdulist = result.master_bias.open()
    assert len(hdulist) == 1

    hdu = hdulist[0]
    assert hdu.shape == (4112, 4096)

    data = hdu.data
    mlevel = 0.0

    block = data[1980:2020, 1980:2020]
    mblock = block.mean()
    sblock = block.std()

    assert abs(mblock - mlevel) < 5 * sblock

    # In the end, remove the files
    for f in fs:
        os.remove(f.name)
