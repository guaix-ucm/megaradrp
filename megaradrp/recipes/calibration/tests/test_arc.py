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

'''Tests for the arc mode recipe module.'''

import os

import pytest

from numina.tests.testcache import download_cache

from numina.user.cli import main
from numina.core import init_drp_system, import_object
from numina.core import ObservationResult
from numina.core import DataFrame
from megaradrp.recipes import BiasRecipe


BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():
    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


@pytest.mark.remote
def test_mode_arc_calibration_set0(numinatpldir):

    run_recipe()
