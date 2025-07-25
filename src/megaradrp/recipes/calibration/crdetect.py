#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.array import combine
from numina.core import Result, Parameter

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import CRMasks
from megaradrp.processing.combine import generate_crmask
import megaradrp.requirements as reqs


class Recipe(MegaraBaseRecipe):
    """Recipe for generation of cosmic ray mask"""

    # Requirements
    method = Parameter(
        'maskscr',
        description='CR masks generation method',
        choices=['maskscr']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for masks generation',
        optional=True
    )
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()

    # Results
    crmasks = Result(CRMasks)

    def run(self, rinput):
        """Execute the recipe."""
        self.logger.info('start MegaraCrDetection recipe')

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        fmethod = getattr(combine, rinput.method)

        img = generate_crmask(
            rinput, flow1,
            method=fmethod,
            method_kwargs=rinput.method_kwargs
        )

        # Update header
        hdr = img[0].header
        self.set_base_headers(hdr)

        result = self.create_result(crmasks=img)
        self.logger.info('end MegaraCrDetection recipe')
        return result
