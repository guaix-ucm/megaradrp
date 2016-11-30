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

"""Acquisition with LCB"""


import numpy as np
from scipy.spatial import KDTree

from numina.core import Product

from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame


class AcquireLCBRecipe(ImageRecipe):
    """Process LCB images."""

    reduced = Product(ProcessedFrame)
    final = Product(ProcessedRSS)

    def run(self, rinput):

        self.logger.info('starting AC LCB reduction')

        reduced2d, reduced1d = super(AcquireLCBRecipe, self).base_run(rinput)
        # rssdata = rss_data[0].data

        self.save_intermediate_img(reduced2d, 'reduced2d.fits')
        self.save_intermediate_img(reduced1d, 'reduced1d.fits')

        do_sky_subtraction = True
        if do_sky_subtraction:
            self.logger.info('start sky subtraction')
            final, origin, sky = self.run_sky_subtraction(reduced1d)
            self.logger.info('end sky subtraction')
        else:
            final =  reduced1d
            origin = final
            sky = final

        fiberconf = self.datamodel.get_fiberconf(final)
        self.logger.debug("LCB configuration is %s", fiberconf.conf_id)
        rssdata = final[0].data
        cut1 = 1000
        cut2 = 3000

        points = [(0, 0)] # Center of fiber 313
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

        npoints = 19
        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        dis_p, idx_p = kdtree.query(points, k=npoints)

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

            coords = np.asarray(coords)
            flux_per_cell = rssdata[colids, cut1:cut2].mean(axis=1)
            flux_per_cell_total = flux_per_cell.sum()
            flux_per_cell_norm = flux_per_cell / flux_per_cell_total
            # centroid
            scf = coords.T * flux_per_cell_norm
            centroid = scf.sum(axis=1)
            self.logger.info('centroid: %s', centroid)
            # central coords
            c_coords = coords - centroid
            scf0 = scf - centroid[:, np.newaxis] * flux_per_cell_norm
            mc2 = np.dot(scf0, c_coords)
            self.logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0,0], mc2[1,1], mc2[0,1])

        return self.create_result(final=final, reduced=reduced2d)
