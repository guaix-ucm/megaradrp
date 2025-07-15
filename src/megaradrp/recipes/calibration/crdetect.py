#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.core import Result, Parameter

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import CRMask


class Recipe(MegaraBaseRecipe):
    """Recipe for cosmic ray detection"""
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )
    crmask = Result(CRMask)

    def run(self, rinput):
        """Execute the recipe."""
        self.logger.info('start cr recipe')

        result = self.create_result(crmask=None)
        self.logger.info('end ce recipe')
        return result
