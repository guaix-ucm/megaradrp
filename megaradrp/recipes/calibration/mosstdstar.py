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

"""MOS Standard Star Image Recipe for Megara"""


from scipy.interpolate import interp1d
import astropy.io.fits as fits

from numina.core import Product
from numina.core.requirements import Requirement

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame, ProcessedSpectrum
from megaradrp.types import ReferenceSpectrumTable, ReferenceExtinctionTable
from megaradrp.types import MasterSensitivity


class MOSStandardRecipe(ImageRecipe):
    """Process MOS Standard Star Recipe.

    This recipe processes a set of images
    obtained in **MOS Stardard Star image** mode and returns
    the total flux of the star.

    See Also
    --------
    megaradrp.recipes.calibration.lcbstdstar.LCBStandardRecipe

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

    The flux of the star is computed by adding the 7 fibers corresponding to the bundle
    containing the star and returned as `star_spectrum`.

    """
    position = Requirement(list, "Position of the reference object", default=(0, 0))
    # nrings = 1
    reference_spectrum = Requirement(ReferenceSpectrumTable, "Spectrum of reference star")
    reference_extinction = Requirement(ReferenceExtinctionTable, "Reference extinction")

    reduced_image = Product(ProcessedFrame)
    final_rss = Product(ProcessedRSS)
    reduced_rss = Product(ProcessedRSS)
    sky_rss = Product(ProcessedRSS)
    star_spectrum = Product(ProcessedSpectrum)
    master_sensitivity = Product(MasterSensitivity)

    def run(self, rinput):

        self.logger.info('starting MOSStandardRecipe reduction')

        reduced2d, rss_data = super(MOSStandardRecipe, self).base_run(rinput)

        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data)
        self.logger.info('end sky subtraction')

        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        # 1 + 6 * Sum_i=0^n =  1 + 3 * n * (n +1)
        # In MOS, only 1 ring around central point
        self.logger.debug('adding %d nrings', 1)
        npoints = 7
        self.logger.debug('adding %d fibers', npoints)

        spectra_pack = self.extract_stars(final, rinput.position, npoints)
        pack = spectra_pack[0]
        # FIXME: include cover1 and cover2
        spectrum, cover1, cover2 = pack
        star_spectrum = fits.PrimaryHDU(spectrum, header=final[0].header)

        star_interp = interp1d(rinput.reference_spectrum[:, 0], rinput.reference_spectrum[:, 1])
        extinc_interp = interp1d(rinput.reference_extinction[:, 0],
                                 rinput.reference_extinction[:, 1])

        sens = self.generate_sensitivity(final, spectrum, star_interp, extinc_interp, cover1, cover2)
        self.logger.info('end MOSStandardRecipe reduction')

        return self.create_result(
            reduced_image=reduced2d,
            final_rss=final,
            reduced_rss=origin,
            sky_rss=sky,
            star_spectrum=star_spectrum,
            master_sensitivity=sens
        )
