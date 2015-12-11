import os

import pytest

from numina.tests.testcache import download_cache

from numina.core import import_object
from numina.core.pipeline import DrpSystem
from numina.core import ObservationResult
from numina.core import DataFrame
from megaradrp.recipes.calibration.bias import BiasRecipe
from megaradrp.loader import load_drp

@pytest.mark.remote
def test_bias(drpmocker):

    BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/BIAS/%s'
    images = ['e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e']

    fs = [download_cache(BASE_URL % i) for i in images]

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=f.name) for f in fs]

    drpmocker.add_drp('MEGARA', load_drp)

    # Here we could directly import the required pipeline,
    # but the idea is to test all the process
    insdrp = DrpSystem().query_by_name(ob.instrument)
    pipeline = insdrp.pipelines.get('default')
    recipe_fqn = pipeline.recipes.get(ob.mode)
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is BiasRecipe

    # TODO: these should be created by a build_recipe_input method
    recipe = BiasRecipe()
    ri = recipe.create_input(obresult=ob)

    result = recipe.run(ri)
    # assert result.qc >= QC.UNKNOWN

    # Checks on the image
    hdulist = result.biasframe.open()
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

if __name__ == "__main__":
    test_bias()