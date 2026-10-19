#
# Copyright 2011-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import types

import astropy.io.fits as fits
import matplotlib
import numpy as np
import pytest
from numina.core import DataFrame

from megaradrp.simulation.actions import simulate_flat
from megaradrp.testing.common import DSHAPE, make_detector, make_obresult, raw_header
from megaradrp.testing.create_ob import generate_bias

pytest_plugins = ["numina.testing.plugins"]

# The recipes draw plots, use a backend without display
matplotlib.use("Agg")

# The simulated detector has a row of dead pixels
DEAD_ROW = 5
DEAD_COLS = 170
# Seed of the noise of the simulated data
SEED = 20261016


@pytest.fixture(scope="session")
def simulated_flats(tmp_path_factory):
    """Simulated flats and master bias of a detector with a row of dead pixels

    The dead pixels are DEAD_COLS pixels of the row DEAD_ROW of the
    quantum efficiency. The flats have a quantum efficiency of 0.8, five of them with a
    source of 5000 counts (low) and five with 40000 counts (high). The
    simulation is slow, the files are shared by the tests, that must not
    modify them.

    The detector draws the noise with the global generator of numpy, it
    is seeded so that the data are always the same, and its state is
    restored at the end.
    """
    basedir = tmp_path_factory.mktemp("simulated")
    state = np.random.get_state()
    np.random.seed(SEED)
    qe = 0.8 * np.ones(DSHAPE)
    qe[DEAD_ROW, 0:DEAD_COLS] = 0.0
    detector = make_detector(qe=qe)

    header = raw_header()

    def write_flats(prefix, source, number=5):
        names = [basedir / f"{prefix}_{idx}.fits" for idx in range(number)]
        for name in names:
            fits.writeto(name, simulate_flat(detector, exposure=1.0, source=source), header=header)
        return names

    low = write_flats("flat_low", 5000.0)
    high = write_flats("flat_high", 40000.0)
    master_bias = basedir / "master_bias.fits"
    generate_bias(detector, 5, basedir).master_bias.frame.writeto(master_bias)
    np.random.set_state(state)
    return types.SimpleNamespace(
        header=header, low=low, high=high, master_bias=master_bias, dead_row=DEAD_ROW, dead_cols=DEAD_COLS
    )


@pytest.fixture(scope="session")
def master_bpm(simulated_flats):
    """The mask of bad pixels of the simulated flats, and the name of its file"""
    from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe

    ob = make_obresult("MegaraBiasImage", simulated_flats.low + simulated_flats.high, simulated_flats.header)
    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(filename=str(simulated_flats.master_bias)))
    result = recipe.run(ri)
    filename = simulated_flats.master_bias.with_name("master_bpm.fits")
    result.master_bpm.frame.writeto(filename)
    return filename
