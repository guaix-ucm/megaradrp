#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''Tests for the calibration module.'''

import os

import pytest

from numina.tests.testcache import download_cache

from numina.user import main
from numina.core import init_drp_system, import_object
from numina.core import ObservationResult
from numina.core import DataFrame
from megaradrp.recipes import BiasRecipe


BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():
    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


def test_recipe1():

    drps = init_drp_system()
    instrument = drps.get('MEGARA')
    pipeline = instrument.pipelines.get('default')
    recipe_fqn = pipeline.recipes.get('bias_image')
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is BiasRecipe


@pytest.mark.remote
def test_recipe2():

    drps = init_drp_system()

    BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/BIAS/%s'
    images = ['e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e',
              'e99d2937d2c29a27c0ba4eebfcf7918e']

    fs = [download_cache(BASE_URL % i) for i in images]

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=f.name) for f in fs]

    instrument = drps.get(ob.instrument)
    pipeline = instrument.pipelines.get('default')
    recipe_fqn = pipeline.recipes.get(ob.mode)
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is BiasRecipe

    # TODO: these should be created by a RecipeInputBuilder
    recipe = BiasRecipe()
    RR = BiasRecipe.RecipeRequirements  # @UndefinedVariable
    ri = RR(obresult=ob)

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


@pytest.mark.remote
def test_mode_bias_set0(numinatpldir):

    run_recipe()


@pytest.mark.remote
def test_mode_bias_set1(numinatpldir):

    run_recipe()


@pytest.mark.remote
def test_mode_fiber_flat_set0(numinatpldir):

    run_recipe()
