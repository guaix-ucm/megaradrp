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


from numina.core import Product
from numina.core.requirements import Requirement

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame, ProcessedSpectrum


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

    reduced_image = Product(ProcessedFrame)
    final_rss = Product(ProcessedRSS)
    reduced_rss = Product(ProcessedRSS)
    sky = Product(ProcessedRSS)
    star_spectrum = Product(ProcessedSpectrum)

    def run(self, rinput):

        self.logger.info('starting MOSStandardRecipe reduction')

        reduced2d, rss_data = super(MOSStandardRecipe, self).base_run(rinput)

        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data)
        self.logger.info('end sky subtraction')

        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring

        # In MOS, only 1 ring around central point
        npoints = 7
        spectrum = self.extract_stars(final, rinput.position, npoints)

        self.logger.info('end MOSStandardRecipe reduction')

        return self.create_result(
            reduced_image=reduced2d,
            final_rss=final,
            reduced_rss=origin,
            sky=sky,
            star_spectrum=spectrum
        )

    def extract_stars(self, final, position, npoints):
        from scipy.spatial import KDTree

        self.logger.info('extracting star')

        fiberconf = self.datamodel.get_fiberconf(final)
        self.logger.debug("MOS configuration is %s", fiberconf.conf_id)
        rssdata = final[0].data

        points = [position]
        fibers = fiberconf.conected_fibers(valid_only=True)
        grid_coords = []
        for fiber in fibers:
            grid_coords.append((fiber.x, fiber.y))
        # setup kdtree for searching
        kdtree = KDTree(grid_coords)

        # Other posibility is
        # query using radius instead
        # radius = 1.2
        # kdtree.query_ball_point(points, k=7, r=radius)

        dis_p, idx_p = kdtree.query(points, k=npoints)
        totals = []
        self.logger.info('Using %d nearest fibers', npoints)
        for diss, idxs, point in zip(dis_p, idx_p, points):
            # For each point
            self.logger.info('For point %s', point)
            colids = []
            coords = []
            for dis, idx in zip(diss, idxs):
                fiber = fibers[idx]
                colids.append(fiber.fibid - 1)
                coords.append((fiber.x, fiber.y))

            flux_total = rssdata[colids].mean(axis=0)
            totals.append(flux_total)
        return totals
