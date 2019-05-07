#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Tests for the calibration module."""

import pytest

from numina.user.cli import main
import numina.drps
from megaradrp.recipes.calibration.bias import BiasRecipe
from megaradrp.loader import load_drp

BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():
    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


def test_recipe1(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)

    insdrp = numina.drps.get_system_drps().query_by_name('MEGARA')
    pipeline = insdrp.pipelines.get('default')
    recipe = pipeline.get_recipe_object('MegaraBiasImage')

    assert isinstance(recipe, BiasRecipe)


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_bias_set0(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)

    run_recipe()


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_bias_set1(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)

    run_recipe()


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_trace_map_set0(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)
    run_recipe()


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_fiber_flat_set0(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)
    run_recipe()
