#
# Copyright 2011-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import logging

from numina.core import BaseRecipe
from numina.types.qc import QC
from numina.core.requirements import ObservationResultRequirement

import megaradrp.core.correctors as cor
from megaradrp.datamodel import MegaraDataModel


class MegaraBaseRecipe(BaseRecipe):
    """Base class for all MEGARA Recipes

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
        """Validate the input of the recipe"""

        import numina.types.multitype
        # print('ATTRS', recipe_input.attrs())
        # print('STORED', recipe_input.stored())
        # print('tag_names', recipe_input.tag_names())

        # Find reference tags
        ref_tags = {}
        for key, val in recipe_input.attrs().items():
            if key == 'obresult':
                ref_tags = val.tags
                break

        # super(MegaraBaseRecipe, self).validate_input(recipe_input)
        # check all the rest against reference tags
        stored = recipe_input.stored()
        attrs = recipe_input.attrs()
        val_results = []
        for key, val in attrs.items():
            if val is None:
                continue
            # If we evaluate query_expr with ref_tags and the tags from the object
            # The result must be true
            req = stored[key]
            rtype = req.type
            # TODO: develop a method to select images that are valid
            # for different filters or insmodes; perhaps a ANY keyword
            tags = rtype.extract_tags(val)
            # FIXME: this should be handled by the type, not with a special case
            if isinstance(rtype, numina.types.multitype.MultiType):
                # get query from first node
                query_expr = rtype.node_type[0].query_expr
            else:
                query_expr = rtype.query_expr

            q2 = query_expr.fill_placeholders(**ref_tags)
            self.logger.debug('type %s with tags %s, expr %s', rtype, tags, q2)
            is_valid = q2.eval(**tags)
            if not is_valid:
                val_results.append((key, q2, tags, ref_tags))
        msg = 'invalid {} with expression {} and tags {} and obs_tags {}'
        for key, q2, tags, ref_tags in val_results:
            self.logger.error(msg.format(key, q2, tags, ref_tags))
        if val_results:
            raise ValueError('Validation error', val_results)

    def run_qc(self, recipe_input, recipe_result):
        """Run Quality Control checks."""
        recipe_result.qc = QC.GOOD
        return recipe_result

    def types_getter(self):
        from megaradrp.ntypes import MasterBias, MasterDark, MasterBPM, MasterSlitFlat
        from megaradrp.ntypes import DiffuseLightCorrection
        imgtypes = [MasterBPM, MasterBias, MasterDark,
                    MasterSlitFlat, DiffuseLightCorrection]
        getters = [cor.get_corrector_bpm, cor.get_corrector_bias,
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
