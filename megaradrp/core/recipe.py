#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

from numina.core import BaseRecipe
from numina.core import DataFrame
from numina.types.qc import QC
from numina.core.requirements import ObservationResultRequirement
import numina.util.flow as flowmod


import megaradrp.core.correctors as cor
from megaradrp.datamodel import MegaraDataModel


class MegaraBaseRecipe(BaseRecipe):
    """Base clase for all MEGARA Recipes


    Parameters
    ----------
    intermediate_results : bool, optional
                           If True, save intermediate results of the Recipe


    Attributes
    ----------

    obresult : ObservationResult, requirement
    logger :
         recipe logger

    datamodel : MegaraDataModel

    """

    obresult = ObservationResultRequirement()
    logger = logging.getLogger('numina.recipes.megara')
    datamodel = MegaraDataModel()

    def validate_input(self, recipe_input):
        """Method to customize recipe input validation.

        See Also
        --------
        numina.core.validator.validate

        """
        self.logger.info('start validating input')
        super(MegaraBaseRecipe, self).validate_input(recipe_input)
        self.logger.info('end validating input')

    def run_qc(self, recipe_input, recipe_result):
        """Run Quality Control checks."""
        recipe_result.qc = QC.GOOD
        return recipe_result

    def types_getter(self):
        from megaradrp.types import MasterBias, MasterDark, MasterBPM, MasterSlitFlat
        imgtypes = [MasterBPM, MasterBias, MasterDark, MasterSlitFlat]
        getters = [cor.get_corrector_bpm, cor.get_corrector_bias,
                   [cor.get_corrector_dark, cor.get_corrector_gain],
                   cor.get_corrector_slit_flat
                   ]
        return imgtypes, getters

    def get_filters(self, imgtypes, getters):
        import collections
        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            self.logger.debug('get_filters, %s  %s', rtype, getter)
            if rtype is None:
                # Unconditional
                if isinstance(getter, collections.Iterable):
                    used_getters.extend(getter)
                else:
                    used_getters.append(getter)
            else:
                # Search
                for key, val in self.RecipeInput.stored().items():
                    if isinstance(val.type, rtype):
                        if isinstance(getter, collections.Iterable):
                            used_getters.extend(getter)
                        else:
                            used_getters.append(getter)
                        break
                else:
                    pass
        return used_getters

    def init_filters_generic(self, rinput, getters_seq, ins):

        meta = self.gather_info(rinput)
        self.logger.debug('obresult info')
        for entry in meta['obresult']:
            self.logger.debug('frame info is %s', entry)

        reduction_flows = []
        for getters in getters_seq:
            correctors = [getter(rinput, meta, ins, self.datamodel) for getter in getters]
            reduction_flow = flowmod.SerialFlow(correctors)
            reduction_flows.append(reduction_flow)

        return reduction_flows

    def init_filters(self, rinput, ins):
        imgtypes1 = [None]
        getters1 = [[cor.get_corrector_overscan, cor.get_corrector_trimming]]
        getters_f1 = self.get_filters(imgtypes1, getters1)
        imgtypes2, getters2 = self.types_getter()
        getters_f2 = self.get_filters(imgtypes2, getters2)

        getters_seq = [getters_f1, getters_f2]
        return self.init_filters_generic(rinput, getters_seq, ins)

    def gather_info(self, recipeinput):
        klass = recipeinput.__class__
        metadata = {}
        for key in klass.stored():
            val = getattr(recipeinput, key)
            if isinstance(val, DataFrame):
                metadata[key] = self.datamodel.gather_info_dframe(val)
            # elif isinstance(val, ObservationResult):
            elif hasattr(val, 'frames'):
                metadata[key] = self.datamodel.gather_info_oresult(val)
            else:
                pass
        return metadata
