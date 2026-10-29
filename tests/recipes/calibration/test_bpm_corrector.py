#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""The correction of bad pixels in the basic processing"""

import astropy.io.fits as fits
import numpy as np
from numina.array import combine
from numina.core import DataFrame
from numina.core.requirements import ObservationResultRequirement
from numina.processing.combine import basic_processing_with_combination

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.testing.common import make_obresult


class DerivedRecipe(MegaraBaseRecipe):
    """Combine the images with the basic processing, keep the result in reduced"""

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_bpm = MasterBPMRequirement()

    def run(self, rinput):
        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        self.reduced = basic_processing_with_combination(rinput, flow1, method=combine.median)
        return self.create_result()


def test_bpm_corrector(simulated_flats, master_bpm):
    ob = make_obresult("MegaraBiasImage", simulated_flats.low, simulated_flats.header)
    recipe = DerivedRecipe()
    ri = recipe.create_input(
        obresult=ob,
        master_bias=DataFrame(filename=str(simulated_flats.master_bias)),
        master_bpm=DataFrame(filename=str(master_bpm)),
    )
    recipe.run(ri)

    data = recipe.reduced[0].data
    bad = fits.getdata(master_bpm) > 0
    # the correction is recorded, and the good pixels have the level of
    # the flats (5000 counts with a quantum efficiency of 0.8)
    assert recipe.reduced[0].header["NUM-BPM"]
    assert np.allclose(np.median(data[~bad]), 4000, rtol=0.01)
