#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Calibration Recipes for Megara"""


from numina.core import Result, Requirement
from numina.types.datatype import ListOfType
from numina.core.requirements import ObservationResultRequirement

import megaradrp.ntypes as typs
from megaradrp.core.recipe import MegaraBaseRecipe


class Recipe(MegaraBaseRecipe):
    """Process Sensitivity Star Recipe.

    This recipe processes a set of images
    processed by the recipes of
    **Standard star with the FIBER MOS** or
    **Standard star with the LCB IFU**
    and returns the sensitivity correction required
    for flux calibration.

    See Also
    --------
    megaradrp.recipes.combined.extinctionstar.Recipe

    """

    # Requirements
    obresult = ObservationResultRequirement()

    master_extinction = Requirement(typs.ReferenceExtinctionTable,
                                    "Atmospheric extinction")
    reference_spectra = Requirement(ListOfType(typs.ReferenceSpectrumTable),
                                    "Reference spectra of Std stars")

    master_sensitivity = Result(typs.MasterSensitivity)

    def run(self, rinput):

        self.logger.info('starting SensivityStarRecipe reduction')

        result = super(Recipe,self).run(rinput)

        return self.create_result(final=result[0], target=result[1], sky=result[2])
