#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

"""LCB Standard Star Image Recipe for Megara"""


import uuid

import numpy
from scipy.interpolate import interp1d
import astropy.io.fits as fits

from numina.core import Product
from numina.core.requirements import Requirement

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame, ProcessedSpectrum
from megaradrp.types import ReferenceSpectrumTable, ReferenceExtinctionTable
from megaradrp.types import MasterSensitivity


class LCBStandardRecipe(ImageRecipe):
    """Process LCB Standard Star Recipe.

    This recipe processes a set of images
    obtained in **LCB Stardard Star image** mode and returns
    the total flux of the star.

    See Also
    --------
    megaradrp.recipes.calibration.mosstdstar.MOSStandardRecipe

    Notes
    -----
    Images provided by `obresult` are trimmed and corrected
    from overscan, bad pixel mask (if `master_bpm` is not None),
    bias, dark current (if `master_dark` is not None) and
    slit-flat (if `master_slitflat` is not None).

    Images thus corrected are then stacked using the median.
    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_traces` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    The sky is subtracted by combining the the fibers marked as `SKY`
    in the fibers configuration. The RSS with sky subtracted is returned ini the
    field `final_rss` of the recipe result.

    The flux of the star is computed by adding summing the fibers in `nrings` around
    the central spaxel containing the star and returned as `star_spectrum`.

    """

    position = Requirement(list, "Position of the reference object", default=(0, 0))
    nrings = Requirement(int, "Number of rings to extract the star", default=3)
    reference_spectrum = Requirement(ReferenceSpectrumTable, "Spectrum of reference star")
    reference_extinction = Requirement(ReferenceExtinctionTable, "Reference extinction")

    reduced_image = Product(ProcessedFrame)
    final_rss = Product(ProcessedRSS)
    reduced_rss = Product(ProcessedRSS)
    sky_rss = Product(ProcessedRSS)
    star_spectrum = Product(ProcessedSpectrum)
    master_sensitivity = Product(MasterSensitivity)

    def run(self, rinput):

        self.logger.info('starting LCBStandardRecipe reduction')

        reduced2d, rss_data = super(LCBStandardRecipe, self).base_run(rinput)

        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data)
        self.logger.info('end sky subtraction')

        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        # 1 + 6 * Sum_i=0^n =  1 + 3 * n * (n +1)
        # Using three rings around central point
        self.logger.debug('adding %d nrings', rinput.nrings)
        npoints = 1 + 3 * rinput.nrings * (rinput.nrings +1)
        self.logger.debug('adding %d fibers', npoints)

        spectrum = self.extract_stars(final, rinput.position, npoints)
        star_spectrum = fits.PrimaryHDU(spectrum[0], header=final[0].header)
        star_interp = interp1d(rinput.reference_spectrum[:,0], rinput.reference_spectrum[:,1])
        extinc_interp = interp1d(rinput.reference_extinction[:, 0],
                               rinput.reference_extinction[:, 1])

        sens = self.generate_sensitivity(final, spectrum, star_interp, extinc_interp)

        return self.create_result(
            reduced_image=reduced2d,
            final_rss=final,
            reduced_rss=origin,
            sky_rss=sky,
            star_spectrum=star_spectrum,
            master_sensitivity=sens
        )

    def generate_sensitivity(self, final, spectrum, star_interp, extinc_interp):

        crpix, wlr0, delt = self.read_wcs(final[0].header)
        wavelen = wlr0 + delt * (numpy.arange(final[0].shape[1]) - crpix)

        airmass = final[0].header['AIRMASS']
        exptime = final[0].header['EXPTIME']

        response_0 = spectrum[0] / exptime
        valid = response_0 > 0
        # In magAB
        # f(Jy) = 3631 * 10^-0.4 mAB
        response_1 = 3631 * numpy.power(10.0, - 0.4 * (star_interp(wavelen) + extinc_interp(wavelen) * airmass))

        # I'm going to filter invalid values anyway
        with numpy.errstate(invalid='ignore', divide='ignore'):
            ratio = response_1 / response_0
            response_2 = numpy.where(valid, ratio, 1.0)

        sens = fits.PrimaryHDU(response_2, header=final[0].header)
        sens.header['uuid'] = str(uuid.uuid1())
        sens.header['tunit'] = ('Jy', "Final units")
        return sens
