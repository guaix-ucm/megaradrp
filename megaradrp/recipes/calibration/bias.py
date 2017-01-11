#
# Copyright 2014-2017 Universidad Complutense de Madrid
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


from numina.core import Product
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterBias
from megaradrp.requirements import MasterBPMRequirement


class BiasRecipe(MegaraBaseRecipe):
    """Process BIAS images and create a MASTER_BIAS product.

    This recipe process a set of bias images obtained in
    *Bias Image* mode and returns a combined product image,
    trimmed to the physical size of the detector.

    Notes
    -----
    Images are corrected from overscan and trimmed to the physical size of the detector.
    Then, they corrected from Bad Pixel Mask, if the BPM is available,
    Finally, images are stacked using the median.

    See Also
    --------
    megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe: recipe to generate MasterBPM


    """

    master_bpm = MasterBPMRequirement()
    master_bias = Product(MasterBias)

    def run(self, rinput):
        """Execute the recipe.

        Parameters
        ----------

        rinput : BiasRecipe.RecipeInput

        Returns
        -------
        BiasRecipe.RecipeResult

        """
        self.logger.info('start bias recipe')
        flow = self.init_filters(rinput, rinput.obresult.configuration)
        errors  = False
        if not errors:
            self.logger.info('not computing errors')
        hdulist = basic_processing_with_combination(rinput, flow, method=combine.median, errors=errors)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        result = self.create_result(master_bias=hdulist)
        self.logger.info('end bias recipe')
        return result
