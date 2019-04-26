#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Dark current calibration Recipe for Megara"""

from numina.core import Result
from numina.array import combine

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.types import MasterDark
from megaradrp.processing.combine import basic_processing_with_combination


class DarkRecipe(MegaraBaseRecipe):

    """Process DARK images and provide MASTER_DARK."""

    master_bias = MasterBiasRequirement()

    master_dark = Result(MasterDark)

    def run(self, rinput):

        flow = self.init_filters(rinput, rinput.obresult.configuration)
        hdulist = basic_processing_with_combination(rinput, flow, method=combine.median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        result = self.create_result(master_dark=hdulist)
        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(DarkRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MasterDark', 'Product type')
        return hdr
