#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

"""Calibration Recipes for Megara"""

from numina.core import Product, Requirement
from numina.core.types import ListOfType
from numina.core.requirements import ObservationResultRequirement

import megaradrp.types as typs
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

    master_extinction = Requirement(typs.Extinction, "Atmospheric extinction")
    reference_spectra = Requirement(ListOfType(typs.ReferenceSpectrum), "Reference spectra of Std stars")

    master_sensitivity = Product(typs.MasterSensitivity)

    def run(self, rinput):

        self.logger.info('starting SensivityStarRecipe reduction')

        result = super(Recipe,self).run(rinput)

        return self.create_result(final=result[0], target=result[1], sky=result[2])
