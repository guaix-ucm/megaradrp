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

import logging

import numpy as np
from astropy.io import fits
from numina.core import BaseRecipe
from numina.core.dataholders import Product
from numina.core.products import QualityControlProduct
from numina.core.requirements import ObservationResultRequirement
from numina.core import DataFrame, ObservationResult
from numina.core.qc import QC
from numina.flow import SerialFlow

import megaradrp.core.correctors as cor
from megaradrp.processing.datamodel import MegaraDataModel


class MegaraBaseRecipe(BaseRecipe):
    """Base clase for all MEGARA Recipes


    Parameters
    ----------
    intermediate_results : bool, optional
                           If True, save intermediate results of the Recipe


    Attributes
    ----------

    obresult : ObservationResult, requirement

    qc : QualityControl, result, QC.GOOD by default

    logger :
         recipe logger

    datamodel : MegaraDataModel

    """

    obresult = ObservationResultRequirement()
    qc = Product(QualityControlProduct, destination='qc', default=QC.UNKNOWN)
    logger = logging.getLogger('numina.recipes.megara')
    datamodel = MegaraDataModel()

    def validate_input(self, recipe_input):
        """Method to customize recipe input validation.

        See Also
        --------
        numina.core.validator.validate

        """
        self.logger.info('start validating input')
        super(MegaraBaseRecipe, self).validate_input(recipe_input)
        self.logger.info('end validating input')

    def run_qc(self, recipe_input, recipe_result):
        """Run Quality Control checks."""
        recipe_result.qc = QC.GOOD
        return recipe_result

    def _resample_rss_flux(self, rss_old, wcalib, wvpar_dict):
        """

        :param rss_old: rss image
        :param wcalib: ndarray of the coefficients
        :param wvpar_dict: dictionary containing wavelength calibration parameters
        :return:
        """
        import math
        from numpy.polynomial.polynomial import polyval
        from numina.array.interpolation import SteffenInterpolator

        nfibers = rss_old.shape[0]
        nsamples = rss_old.shape[1]

        # print nfibers, nsamples
        # z = [0, nsamples - 1]
        # res = polyval(z, wcalib.T)
        # print res
        # all_delt = (res[:, 1] - res[:, 0]) / nsamples
        # print all_delt.max(), all_delt.min(), np.median(all_delt)
        #
        # delts = all_delt.min()
        # delts = np.median(all_delt)
        # delts = 0.37
        # print 'median of delts', delts
        #
        # # first pixel is
        # wl_min = res[:, 0].min()
        # wl_min = 7140.0 #res[:, 0].min()
        # # last pixel is
        # wl_max = res[:, 1].max()
        # wl_max = 8730.63
        # print 'pixel range', wl_min, wl_max
        #
        # npix = int(math.ceil((wl_max - wl_min) / delts))

        npix = wvpar_dict['npix']
        delts = wvpar_dict['cdelt']
        wl_min = wvpar_dict['crval']
        crpix = wvpar_dict['crpix']

        wl_max = wl_min + (npix - crpix) * delts

        new_x = np.arange(npix)
        new_wl = wl_min + delts * new_x

        old_x_borders = np.arange(-0.5, nsamples)
        old_x_borders += crpix  # following FITS criterium
        old_wl_borders = polyval(old_x_borders, wcalib.T)

        new_borders = self._map_borders(new_wl)

        accum_flux = np.empty((nfibers, nsamples + 1))
        accum_flux[:, 1:] = np.cumsum(rss_old, axis=1)
        accum_flux[:, 0] = 0.0
        rss_resampled = np.zeros((nfibers, npix))

        for idx  in range(nfibers):
            # We need a monotonic interpolator
            # linear would work, we use a cubic interpolator
            interpolator = SteffenInterpolator(old_wl_borders[idx],accum_flux[idx], extrapolate='border')
            fl_borders = interpolator(new_borders)
            rss_resampled[idx] = fl_borders[1:] - fl_borders[:-1]

        return rss_resampled, (wl_min, wl_max, delts)

    def _map_borders(self, wls):
        """Compute borders of pixels for interpolation.

        The border of the pixel is assumed to be midway of the wls
        """
        midpt_wl = 0.5 * (wls[1:] + wls[:-1])
        all_borders = np.zeros((wls.shape[0] + 1,))
        all_borders[1:-1] = midpt_wl
        all_borders[0] = 2 * wls[0] - midpt_wl[0]
        all_borders[-1] = 2 * wls[-1] - midpt_wl[-1]
        return all_borders

    def add_wcs(self, hdr, wlr0, delt, crpix=1.0):
        hdr['CRPIX1'] = crpix
        hdr['CRVAL1'] = wlr0
        hdr['CDELT1'] = delt
        hdr['CTYPE1'] = 'WAVELENGTH'
        hdr['CRPIX2'] = 1
        hdr['CRVAL2'] = 1
        hdr['CDELT2'] = 1
        hdr['CTYPE2'] = 'PIXEL'
        return hdr

    def types_getter(self):
        from megaradrp.types import MasterBias, MasterDark, MasterBPM, MasterSlitFlat
        imgtypes = [None, MasterBPM, MasterBias, MasterDark, MasterSlitFlat]
        getters = [[cor.get_corrector_overscan, cor.get_corrector_trimming],
                   cor.get_corrector_bpm, cor.get_corrector_bias,
                   [cor.get_corrector_dark, cor.get_corrector_gain],
                   cor.get_corrector_slit_flat
                   ]
        return imgtypes, getters

    def get_filters(self):
        import collections
        imgtypes, getters = self.types_getter()
        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            self.logger.debug('get_filters, %s  %s', rtype, getter)
            if rtype is None:
                # Unconditional
                if isinstance(getter, collections.Iterable):
                    used_getters.extend(getter)
                else:
                    used_getters.append(getter)
            else:
                # Search
                for key, val in self.RecipeInput.stored().items():
                    if isinstance(val.type, rtype):
                        if isinstance(getter, collections.Iterable):
                            used_getters.extend(getter)
                        else:
                            used_getters.append(getter)
                        break
                else:
                    pass
        return used_getters

    def init_filters_generic(self, rinput, getters, ins):

        meta = self.gather_info(rinput)
        self.logger.debug('obresult info')
        for entry in meta['obresult']:
            self.logger.debug('frame info is %s', entry)
        correctors = [getter(rinput, meta, ins, self.datamodel) for getter in getters]

        flow = SerialFlow(correctors)

        return flow

    def init_filters(self, rinput, ins):
        getters = self.get_filters()
        return self.init_filters_generic(rinput, getters, ins)

    def gather_info(self, recipeinput):
        klass = recipeinput.__class__
        metadata = {}
        for key in klass.stored():
            val = getattr(recipeinput, key)
            if isinstance(val, DataFrame):
                metadata[key] = self.datamodel.gather_info_dframe(val)
            elif isinstance(val, ObservationResult):
                metadata[key] = self.datamodel.gather_info_oresult(val)
            else:
                pass
        return metadata
