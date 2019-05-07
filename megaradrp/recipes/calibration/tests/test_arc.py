#
# Copyright 2015-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Tests for the arc mode recipe module."""

import pytest

from numina.user.cli import main

from megaradrp.loader import load_drp


BASE_URL = 'http://guaix.fis.ucm.es/~spr/megara_test/'


def run_recipe():

    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])


@pytest.mark.remote_data
@pytest.mark.usefixtures("numinatpldir")
def test_mode_arc_calibration_set0(drpmocker):

    drpmocker.add_drp('MEGARA', load_drp)

    run_recipe()