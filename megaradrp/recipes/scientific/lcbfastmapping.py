#
# Copyright 2011-201t Universidad Complutense de Madrid
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
