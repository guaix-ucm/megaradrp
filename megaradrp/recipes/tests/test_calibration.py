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

from numina.user.cli import main
from numina.core import import_object
from numina.core.pipeline import DrpSystem
from megaradrp.recipes.calibration.BiasRecipe import BiasRecipe
from megaradrp.loader import megara_drp_load, load_cli_storage

BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():
    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


def test_recipe1(drpmocker):

    drpmocker.add_drp('MEGARA', megara_drp_load)

    insdrp = DrpSystem().query_by_name('MEGARA')
    pipeline = insdrp.pipelines.get('default')

    recipe_fqn = pipeline.recipes.get('bias_image')
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is BiasRecipe

@pytest.mark.remote
@pytest.mark.usefixtures("numinatpldir")
def test_mode_bias_set0(drpmocker):

    drpmocker.add_drp('MEGARA', megara_drp_load)

    run_recipe()


@pytest.mark.remote
@pytest.mark.usefixtures("numinatpldir")
def test_mode_bias_set1(drpmocker):
    drpmocker.add_drp('MEGARA', megara_drp_load)

    run_recipe()

@pytest.mark.remote
@pytest.mark.usefixtures("numinatpldir")
def test_mode_trace_map_set0(drpmocker):
    drpmocker.add_drp('MEGARA', megara_drp_load)

    run_recipe()

@pytest.mark.remote
@pytest.mark.usefixtures("numinatpldir")
def test_mode_fiber_flat_set0(drpmocker):
    drpmocker.add_drp('MEGARA', megara_drp_load)

    # FIXME: this was completely crazy to debug
    # if we don't load store/dump functions
    # tracemap is loaded as a its filename
    load_cli_storage()

    run_recipe()