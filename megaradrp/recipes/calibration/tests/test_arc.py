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

"""Tests for the arc mode recipe module."""

import pytest

from numina.user.cli import main

from megaradrp.loader import load_drp


BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():

    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


@pytest.mark.remote
@pytest.mark.usefixtures("numinatpldir")
def test_mode_arc_calibration_set0(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)

    run_recipe()