#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""LCB Fast Mapping Recipe for Megara"""


from numina.core import Product, ObservationResult

from megaradrp.types import ProcessedMultiRSS
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.processing.multirss import generate_multi_rss


class LCBFastMappingRecipe(MegaraBaseRecipe):
    """Process LCB Fast Mapping Recipe.

    This recipe processes a set of images
    obtained in **LCB Fast Mapping image** mode and
    a multi RSS image.

    See Also
    --------
    megaradrp.recipes.scientific.lcb.LCBImageRecipe
    megaradrp.recipes.scientific.mos.MOSImageRecipe

    Notes
    -----
    Images previously obtained in **LCB image** and reduced
    are stacked together in multi RSS format.

    """
    final_multirss = Product(ProcessedMultiRSS)

    def run(self, rinput):
        self.logger.info('start FastMappingRecipe')
        obresult = rinput.obresult
        imgs = [frame.open() for frame in obresult.frames]

        result = generate_multi_rss(imgs)
        self.logger.info('end FastMappingRecipe')

        return self.create_result(final_multirss=result)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        """Quey previous results of LCB image"""
        return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):
        cls.logger.debug('start recipe input builder')
        stareImagesIds = obsres.stareImagesIds
        cls.logger.debug('LCB image IDS %s: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['final_rss'])

        newOR = ObservationResult()
        newOR.frames = stareImages
        newRI = cls.create_input(obresult=newOR)
        cls.logger.debug('end recipe input builder')
        return newRI
