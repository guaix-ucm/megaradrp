#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import astropy.io.fits as fits
import numpy
from numina.core import DataFrame

from megaradrp.recipes.calibration.dark import DarkRecipe
from megaradrp.simulation.actions import simulate_dark_fits
from megaradrp.simulation.factory import MegaraImageFactory
from megaradrp.testing.common import DSHAPE, make_detector, make_obresult, raw_header


def test_dark(tmp_path):
    numpy.random.seed(422992983)

    dark = 3.0  # In 1 hour
    exptime = 3600.0
    detector = make_detector(qe=0.8 * numpy.ones(DSHAPE), dark=dark / exptime, ron=0.001)

    number = 3
    names = []
    for idx, img in enumerate(simulate_dark_fits(MegaraImageFactory(), detector, exposure=exptime, repeat=number)):
        name = tmp_path / f"dark_{idx}.fits"
        img.writeto(name)
        names.append(name)

    header = raw_header()
    master_bias_hdul = fits.HDUList(fits.PrimaryHDU(numpy.zeros(DSHAPE), header=header))

    ob = make_obresult("MegaraDarkImage", names, header)
    recipe = DarkRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(frame=master_bias_hdul))
    result = recipe.run(ri)

    mean_dark_value = result.master_dark.frame[0].data.mean()
    assert numpy.allclose(mean_dark_value, dark, atol=0, rtol=1e-1)
