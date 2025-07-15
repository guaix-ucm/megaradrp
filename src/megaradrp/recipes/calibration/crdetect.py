#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.core import Result, Parameter
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import MasterBias
from megaradrp.requirements import MasterBPMRequirement


class Recipe(MegaraBaseRecipe):
    """Recipe for cosmic ray detection"""
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )

    def run(self, rinput):
        """Execute the recipe."""
        self.logger.info('start cr recipe')

        result = self.create_result()
        self.logger.info('end ce recipe')
        return result
