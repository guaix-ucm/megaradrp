#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

"""Dark current calibration Recipe for Megara"""

from numina.core import Product

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement
from megaradrp.products import MasterDark


class DarkRecipe(MegaraBaseRecipe):

    '''Process DARK images and provide MASTER_DARK. '''

    master_bias = MasterBiasRequirement()

    master_dark = Product(MasterDark)

    def __init__(self):
        super(DarkRecipe, self).__init__(version="0.1.0")


    def run(self, rinput):

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        hdu, data = self.hdu_creation(rinput.obresult, {'biasmap':mbias})

        try:
            hdu[0].data = hdu[0].data/hdu[0].header['EXPTIME']
        except:
            pass

        result = self.create_result(master_dark=hdu)
        return result
