#
# Copyright 2011-2016 Universidad Complutense de Madrid
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


import numpy
from numina.core import Product

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import MasterFiberFlat


class LCBImageRecipe(ImageRecipe):
    """Process LCB images."""

    final = Product(MasterFiberFlat)
    target = Product(MasterFiberFlat)
    sky = Product(MasterFiberFlat)

    def __init__(self):
        super(LCBImageRecipe, self).__init__()

    def run(self, rinput):

        self.logger.info('starting LCB reduction')

        reduced, rss_data = super(LCBImageRecipe,self).run(rinput)

        # Take a look at == []
        indices = []
        wlcalib = []
        for key, val in rinput.wlcalib.contents.items():
            if val.coeff:
                wlcalib.append(val.coeff)
                if len(indices)==0:
                    indices.append(0)
                else:
                    indices.append(indices[-1])
            else:
                indices.append(indices[-1]+1)
        wlcalib_aux = numpy.asarray(wlcalib)
        #final, wcsdata = self.resample_rss_flux(rss, wlcalib_aux, indices)
        final, wcsdata = self.resample_rss_flux(rss_data.data, wlcalib_aux, indices)

        #
        import astropy.io.fits as fits
        rss = fits.PrimaryHDU(data=final)
        #

        return self.create_result(final=rss, target=reduced, sky=reduced)
