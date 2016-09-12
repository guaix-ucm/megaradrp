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

'''Calibration Recipes for Megara'''

import logging
from astropy.io import fits
import numpy as np
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement, Requirement
# from numina.flow import SerialFlow

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.types import MasterFiberFlat, WavelengthCalibration
from megaradrp.types import MasterWeights, TraceMap
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.requirements import MasterDarkRequirement, MasterFiberFlatRequirement
from megaradrp.requirements import MasterSlitFlatRequirement, MasterTwilightRequirement
# from megaradrp.processing.fiberflat import FiberFlatCorrector
# from megaradrp.processing.twilight import TwilightCorrector
# from megaradrp.processing.weights import WeightsCorrector
from megaradrp.core.processing import apextract_tracemap_2


_logger = logging.getLogger('numina.recipes.megara')


class ImageRecipe(MegaraBaseRecipe):
    """Base Image."""

    # Requirements  
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_bpm = MasterBPMRequirement()
    master_slitflat = MasterSlitFlatRequirement()
    wlcalib = Requirement(WavelengthCalibration, 'Wavelength calibration table')
    # master_weights = Requirement(MasterWeights, 'Set of files')
    master_fiberflat = MasterFiberFlatRequirement()
    master_twilight = MasterTwilightRequirement()
    tracemap = Requirement(TraceMap, 'Trace information of the Apertures')

    def __init__(self):
        super(ImageRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):

        parameters = self.get_parameters(rinput)
        reduced = self.bias_process_common(rinput.obresult, parameters)
        rssdata = apextract_tracemap_2(reduced[0].data, rinput.tracemap)

        return reduced, rssdata


    def get_wcallib(self, lambda1, lambda2, fibras, traces, rss, neigh_info, grid):


        # Take a look at == []
        indices = []
        wlcalib = []
        for elem in traces:
            if elem['aperture']['function']['coefficients']:
                wlcalib.append(elem['aperture']['function']['coefficients'])
                if len(indices)==0:
                    indices.append(0)
                else:
                    indices.append(indices[-1])
            else:
                indices.append(indices[-1]+1)

        wlcalib_aux = np.asarray(wlcalib)
        final, wcsdata = self.resample_rss_flux(rss, wlcalib_aux, indices)

        hdu_f = fits.PrimaryHDU(final)
        hdu_f.writeto('resample_rss.fits', clobber=True)

        fibras.sort()

        suma = np.sum(final[fibras.astype(int),lambda1:lambda2],axis=1)
        sky = np.min(suma)
        sumaparcial = suma - sky

        neigh_info = neigh_info[np.argsort(neigh_info[:, 0])]

        centroid_x = np.multiply(sumaparcial,neigh_info[:,2])
        centroid_x = np.sum(centroid_x, axis=0)

        centroid_y = np.multiply(sumaparcial,neigh_info[:,3])
        centroid_y = np.sum(centroid_y, axis=0)

        sumatotal = np.sum(sumaparcial, axis=0)
        _logger.info( "total sum: %s", sumatotal)

        second_order = []
        aux = np.sum(np.multiply(suma,(neigh_info[:,2] - np.mean(neigh_info[:,2]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        _logger.info("Second order momentum X: %s", second_order[0])

        aux = np.sum(np.multiply(suma,(neigh_info[:,3] - np.mean(neigh_info[:,3]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        _logger.info("Second order momentum Y: %s", second_order[1])

        aux = np.multiply(neigh_info[:,3] - np.mean(neigh_info[:,3]),neigh_info[:,2] - np.mean(neigh_info[:,2]))
        aux = np.sum(np.multiply(aux,suma))
        cov = np.divide(aux ,np.sum(suma, axis=0))
        _logger.info("Cov X,Y: %s", cov)

        centroid_x = np.divide(centroid_x, sumatotal)
        _logger.info( "centroid_x: %s", centroid_x)

        centroid_y = np.divide(centroid_y, sumatotal)
        _logger.info("centroid_y: %s", centroid_y)

        centroid = [centroid_x, centroid_y]

        peak = np.sum(final[grid.get_fiber(centroid),lambda1:lambda2],axis=0)

        return centroid, sky, peak, second_order, cov


    def generate_solution(self, points, centroid, sky, fiber, peaks, second_order, cova):
        result = []
        for cont, value in enumerate(points):
            lista = (value[0], value[1], centroid[cont][0],centroid[cont][1], sky[cont], fiber[cont], peaks[cont], second_order[cont][0], second_order[cont][1], cova[cont])
            result.append(lista)
        return np.array(result, dtype=[('x_point','float'),('y_point','float'),('x_centroid','float'),('y_centroid','float'), ('sky','float'),('fiber','int'),('peak','float'),('x_second_order','float'), ('y_second_order','float'), ('covariance','float') ])

    def generateJSON(self, points, centroid, sky, fiber, peaks, second_order, cova):
        '''
        '''

        _logger.info('start JSON generation')

        result = []
        for cont, value in enumerate(points):
            obj = {
                'points': value,
                'centroid': centroid[cont],
                'sky':sky[cont],
                'fiber': fiber[cont],
                'peak': peaks[cont],
                'second_order': second_order[cont],
                'covariance': cova[cont]
            }
            result.append(obj)

        _logger.info('end JSON generation')

        return result