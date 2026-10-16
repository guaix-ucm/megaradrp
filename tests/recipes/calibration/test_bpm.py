#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Tests for the bpm mode recipe module."""

import astropy.io.fits as fits
import numpy as np


def test_bpm(simulated_flats, master_bpm):
    """The mask of the simulated flats, computed by the fixture master_bpm"""
    bpm = fits.getdata(master_bpm)
    ncols = simulated_flats.dead_cols

    # the dead pixels are marked, in one row, and no other pixel
    rows, cols = np.nonzero(bpm)
    assert len(set(rows)) == 1
    assert cols.min() >= 0 and cols.max() < ncols
    assert ncols - 10 <= len(cols) <= ncols
