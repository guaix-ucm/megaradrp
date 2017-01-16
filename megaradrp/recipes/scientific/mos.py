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

"""MOS Image Recipe for Megara"""

from numina.core import Product

from megaradrp.types import ProcessedRSS, ProcessedFrame
from .base import ImageRecipe


class MOSImageRecipe(ImageRecipe):
    """Process MOS images.

    This recipe processes a set of images
    obtained in **MOS image** mode and returns
    the sky subtracted RSS.

    See Also
    --------
    megaradrp.recipes.scientific.lcb.LCBImageRecipe

    Notes
    -----
    Images provided by `obresult` are trimmed and corrected
    from overscan, bad pixel mask (if `master_bpm` is not None),
    bias, dark current (if `master_dark` is not None) and
    slit-flat (if `master_slitflat` is not None).

    Images thus corrected are the stacked using the median.
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

    """

    reduced_image = Product(ProcessedFrame)
    final_rss = Product(ProcessedRSS)
    reduced_rss = Product(ProcessedRSS)
    sky = Product(ProcessedRSS)

    def run(self, rinput):
        self.logger.info('starting MOS reduction')

        reduced2d, rss_data = super(MOSImageRecipe, self).base_run(rinput)

        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data)
        self.logger.info('end sky subtraction')
        self.logger.info('end MOS reduction')

        return self.create_result(
            reduced_image=reduced2d,
            final_rss=final,
            reduced_rss=origin,
            sky=sky
        )
