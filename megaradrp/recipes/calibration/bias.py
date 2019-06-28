#
# Copyright 2014-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.core import Result
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterBias
from megaradrp.requirements import MasterBPMRequirement


class BiasRecipe(MegaraBaseRecipe):
    """Process BIAS images and create a MASTER_BIAS product.

    This recipe process a set of bias images obtained in
    **Bias Image** mode and returns a combined product image,
    trimmed to the physical size of the detector.

    Notes
    -----
    Images are corrected from overscan and trimmed to the physical size of the detector.
    Then, they corrected from Bad Pixel Mask, if the BPM is available,
    Finally, images are stacked using the median.

    See Also
    --------
    megaradrp.types.MasterBias: description of the MasterBias product

    """

    master_bpm = MasterBPMRequirement()
    master_bias = Result(MasterBias)

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

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(BiasRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterBias', 'Product type')
        return hdr
