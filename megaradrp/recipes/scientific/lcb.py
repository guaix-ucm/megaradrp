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


from astropy.io import fits
import numpy
from numina.core import Product

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import MasterFiberFlat

LCB_NFIBERS = 623

# FIXME: hardcoded numbers
vph_thr = {'default': {'LR-I':{'crval': 7140.0,
                              'cdelt': 0.37,
                              'crpix': 1.0,
                              'npix': 4300},
                      },
}


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

        rss_data.data = numpy.fliplr(rss_data.data)
        current_vph = rinput.obresult.tags['vph']
        if current_vph not in vph_thr['default']:
            raise ValueError('grism ' + current_vph + ' is not defined in ' +
                             'vph_thr dictionary')
        wvpar_dict = vph_thr['default'][current_vph]

        wlcalib = []
        for fidx in range(1, LCB_NFIBERS + 1):
            if fidx in rinput.wlcalib.contents:
                wlcalib.append(rinput.wlcalib.contents[fidx].coeff)
            else:
                # FIXME: polynomial degree forced to be 5
                wlcalib.append(numpy.array([0., 1., 0., 0., 0., 0.]))

        wlcalib_aux = numpy.asarray(wlcalib)
        final, wcsdata = self.resample_rss_flux(rss_data.data,
                                                wlcalib_aux, wvpar_dict)

        rss = fits.PrimaryHDU(data=final.astype(numpy.float32),
                              header=rss_data.header)
        self.add_wcs(rss.header, wvpar_dict['crval'], wvpar_dict['cdelt'],
                     wvpar_dict['crpix'])

        return self.create_result(final=rss, target=reduced, sky=reduced)
