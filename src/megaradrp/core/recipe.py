#
# Copyright 2011-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.core import BaseRecipe
from numina.core.recipeinout import RecipeInput, TagsChecker
from numina.types.qc import QC
from numina.core.requirements import ObservationResultRequirement

import megaradrp.core.correctors as cor
from megaradrp.datamodel import MegaraDataModel


class MegaraRecipeInput(RecipeInput):
    """Input of the MEGARA recipes

    The tags of the products of the input must match those of the
    observation result.
    """

    __checkers__ = [TagsChecker()]


class MegaraBaseRecipe(BaseRecipe):
    """Base class for all MEGARA Recipes

    Each recipe logs with the logger of its module.

    Attributes
    ----------
    obresult : ObservationResult, requirement
    datamodel : MegaraDataModel
    """

    obresult = ObservationResultRequirement()
    datamodel = MegaraDataModel()
    RecipeInput = MegaraRecipeInput

    def run_qc(self, recipe_input, recipe_result):
        """Run Quality Control checks."""
        recipe_result.qc = QC.GOOD
        return recipe_result

    def types_getter(self):
        from megaradrp.ntypes import MasterBias, MasterDark, MasterBPM, MasterSlitFlat
        from megaradrp.ntypes import DiffuseLightCorrection

        imgtypes = [MasterBPM, MasterBias, MasterDark, MasterSlitFlat, DiffuseLightCorrection]
        getters = [
            cor.get_corrector_bpm,
            cor.get_corrector_bias,
            [cor.get_corrector_dark, cor.get_corrector_gain],
            cor.get_corrector_slit_flat,
            cor.get_corrector_diffuse_light,
        ]
        return imgtypes, getters

    def create_getters(self):
        imgtypes1 = [None]
        getters1 = [[cor.get_corrector_overscan, cor.get_corrector_trimming]]
        getters_f1 = self.get_filters(imgtypes1, getters1)
        imgtypes2, getters2 = self.types_getter()
        getters_f2 = self.get_filters(imgtypes2, getters2)
        getters_seq = [getters_f1, getters_f2]
        return getters_seq
