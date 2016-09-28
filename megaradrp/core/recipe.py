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

import logging

import numina.array.combine as combine
import numpy as np
from astropy.io import fits
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector, BadPixelCorrector
from numina.flow.processing import DarkCorrector
from numina.core import BaseRecipe
from numina.core.dataholders import Product
from numina.core.products import QualityControlProduct
from numina.core.requirements import ObservationResultRequirement

import megaradrp.core.correctors as cor
from megaradrp.processing.trimover import OverscanCorrector, TrimImage
from megaradrp.processing.slitflat import SlitFlatCorrector
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.fiberflat import FiberFlatCorrector
from megaradrp.processing.datamodel import MegaraDataModel
from megaradrp.processing.twilight import TwilightCorrector


class MegaraBaseRecipe(BaseRecipe):
    """Base clase for all MEGARA Recipes"""

    obresult = ObservationResultRequirement()
    qc = Product(QualityControlProduct, dest='qc')
    logger = logging.getLogger('numina.recipes.megara')

    def __init__(self, version):
        self.__flow = {'ArcCalibrationRecipe': [OverscanCorrector, TrimImage,
                                                BiasCorrector,
                                                BadPixelCorrector,
                                                DarkCorrector],
                       'BadPixelsMaskRecipe': [OverscanCorrector, TrimImage,
                                               BiasCorrector, DarkCorrector],
                       'FiberFlatRecipe': [OverscanCorrector, TrimImage,
                                           BiasCorrector, DarkCorrector,
                                           BadPixelCorrector,
                                           SlitFlatCorrector],
                       'SlitFlatRecipe': [OverscanCorrector, TrimImage,
                                          BiasCorrector, BadPixelCorrector,
                                          DarkCorrector],
                       'TraceMapRecipe': [OverscanCorrector, TrimImage,
                                          BiasCorrector, BadPixelCorrector,
                                          DarkCorrector],
                       'WeightsRecipe': [OverscanCorrector, TrimImage,
                                         BiasCorrector, BadPixelCorrector,
                                         DarkCorrector, SlitFlatCorrector],
                       'TwilightFiberFlatRecipe': [OverscanCorrector,
                                                   TrimImage,
                                                   BiasCorrector,
                                                   BadPixelCorrector,
                                                   DarkCorrector,
                                                   SlitFlatCorrector],
                       'AcquireLCBRecipe': [OverscanCorrector, TrimImage,
                                          BiasCorrector, BadPixelCorrector,
                                          DarkCorrector, SlitFlatCorrector],
                                          # WeightsCorrector,
                                          # FiberFlatCorrector,
                                          # TwilightCorrector],
                       'AcquireMOSRecipe': [OverscanCorrector, TrimImage,
                                          BiasCorrector, BadPixelCorrector,
                                          DarkCorrector, SlitFlatCorrector],
                                          # WeightsCorrector,
                                          # FiberFlatCorrector,
                                          # TwilightCorrector],
                       'PseudoFluxCalibrationRecipe': [OverscanCorrector,
                                                       TrimImage,
                                                       BiasCorrector,
                                                       BadPixelCorrector,
                                                       DarkCorrector,
                                                       ApertureExtractor,
                                                       FiberFlatCorrector],
                       }
        super(MegaraBaseRecipe, self).__init__(version=version)

    def __generate_flow(self, params, confFile):
        import copy
        ff = self.__flow[self.__class__.__name__]
        flow = copy.deepcopy(ff)
        try:
            cont = 0
            while cont < len(flow):
                if issubclass(BiasCorrector, flow[cont]):
                    flow[cont] = (flow[cont](params['biasmap']))
                elif issubclass(BadPixelCorrector, flow[cont]):
                    if 'bpm' in params.keys():
                        flow[cont] = (flow[cont](params['bpm']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(DarkCorrector, flow[cont]):
                    if 'dark' in params.keys():
                        flow[cont] = (flow[cont](params['dark']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(SlitFlatCorrector, flow[cont]):
                    if 'slitflat' in params.keys():
                        flow[cont] = (flow[cont](params['slitflat']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(ApertureExtractor, flow[cont]):
                    if 'traces' in params.keys():
                        flow[cont] = (flow[cont](params['traces']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(FiberFlatCorrector, flow[cont]):
                    if 'fiberflat' in params.keys():
                        flow[cont] = (flow[cont](params['fiberflat']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(TwilightCorrector, flow[cont]):
                    if 'twilight' in params.keys():
                        flow[cont] = (flow[cont](params['twilight']))
                    else:
                        del (flow[cont])
                        cont -= 1
                elif issubclass(TrimImage, flow[cont]) or issubclass(
                        OverscanCorrector, flow[cont]):
                    flow[cont] = (flow[cont](confFile=confFile))
                cont += 1
            basicflow = SerialFlow(flow)

        except Exception as e:
            self.logger.error(e)
            raise (e)
        del flow
        return basicflow

    def bias_process_common(self, obresult, img):

        hdu, data = self.hdu_creation(obresult, img)

        hdr = hdu[0].header
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()
        if 'FILENAME' in hdr.keys():
            del hdr['FILENAME']

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        reduced = fits.HDUList(hdu + [varhdu, num])
        return reduced

    def hdu_creation(self, obresult, params=None):

        if params is None:
            params = {}

        basicflow = self.__generate_flow(params, obresult.configuration.values)
        cdata = []
        headers = []
        hdulist = []
        try:
            for frame in obresult.images:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                cdata.append(hdulist)

            self.logger.info('stacking %d images using median', len(cdata))

            data = combine.median([d[0].data for d in cdata], dtype='float32')
            template_header = cdata[0][0].header

            for header in hdulist:
                if 'PRIMARY' in header.name:
                    headers.append(fits.PrimaryHDU(data[0], header=template_header))
                else:
                    headers.append(header)

        finally:
            for hdulist in cdata:
                hdulist.close()

        return fits.HDUList(headers), data

    def get_parameters(self, rinput):

        parameters = {}

        try:
            if rinput.master_bias:
                with rinput.master_bias.open() as hdul:
                    parameters['biasmap'] = hdul[0].data.copy()
        except:
            pass

        try:
            if rinput.master_bpm:
                with rinput.master_bpm.open() as hdul:
                    parameters['bpm'] = hdul[0].data.copy()
        except:
            pass
        try:
            if rinput.master_dark:
                with rinput.master_dark.open() as hdul:
                    parameters['dark'] = hdul[0].data.copy()
        except:
            pass

        try:
            if rinput.master_slitflat:
                with rinput.master_slitflat.open() as hdul:
                    parameters['slitflat'] = hdul[0].data.copy()
        except:
            pass

        try:
            if rinput.master_fiberflat:
                with rinput.master_fiberflat.open() as hdul:
                    parameters['fiberflat'] = hdul[0].data.copy()
        except:
            pass

        try:
            if rinput.master_twilight:
                with rinput.master_twilight.open() as hdul:
                    parameters['twilight'] = hdul[0].data.copy()
        except:
            pass

        try:
            if rinput.master_weights:
                parameters['weights'] = rinput.master_weights
        except:
            pass

        try:
            if rinput.traces:
                parameters['traces'] = rinput.traces
        except:
            pass

        return parameters

    def get_wlcalib(self, data):

        wlcalib = [elem['aperture']['function']['coefficients'] for elem in
                   data]
        return np.array(wlcalib)

    def resample_rss_flux(self, rss_old, wcalib, indexes=None):
        """

        :param rss_old: rss image
        :param wcalib: ndarray of the coefficients
        :param indexes: is an array to take into account that some wlcalib might not be done
        :return:
        """
        import math
        from numpy.polynomial.polynomial import polyval
        from numina.array.interpolation import SteffenInterpolator

        if indexes is None:
            indexes = []

        nfibers = rss_old.shape[0]
        nsamples = rss_old.shape[1]
        z = [0, nsamples - 1]
        res = polyval(z, wcalib.T)
        all_delt = (res[:, 1] - res[:, 0]) / nsamples

        delts = all_delt.min()

        # first pixel is
        wl_min = res[:, 0].min()
        # last pixel is
        wl_max = res[:, 1].max()

        npix = int(math.ceil((wl_max - wl_min) / delts))
        new_x = np.arange(npix)
        new_wl = wl_min + delts * new_x

        old_x_borders = np.arange(-0.5, nsamples)
        old_wl_borders = polyval(old_x_borders, wcalib.T)

        new_borders = self.map_borders(new_wl)

        accum_flux = np.empty((nfibers, nsamples + 1))
        accum_flux[:, 1:] = np.cumsum(rss_old, axis=1)
        accum_flux[:, 0] = 0.0
        rss_resampled = np.zeros((nfibers, npix))
        for idx in range(nfibers):
            # We need a monotonic interpolator
            # linear would work, we use a cubic interpolator
            if len(indexes)==0:
                interpolator = SteffenInterpolator(old_wl_borders[idx],accum_flux[idx],extrapolate='border')
            else:
                interpolator = SteffenInterpolator(old_wl_borders[idx-indexes[idx]],accum_flux[idx],extrapolate='border')
            fl_borders = interpolator(new_borders)
            rss_resampled[idx] = fl_borders[1:] - fl_borders[:-1]
        return rss_resampled, (wl_min, wl_max, delts)

    def map_borders(self, wls):
        """Compute borders of pixels for interpolation.

        The border of the pixel is assumed to be midway of the wls
        """
        midpt_wl = 0.5 * (wls[1:] + wls[:-1])
        all_borders = np.zeros((wls.shape[0] + 1,))
        all_borders[1:-1] = midpt_wl
        all_borders[0] = 2 * wls[0] - midpt_wl[0]
        all_borders[-1] = 2 * wls[-1] - midpt_wl[-1]
        return all_borders

    def add_wcs(self, hdr, wlr0, delt):
        hdr['CRPIX1'] = 1
        hdr['CRVAL1'] = wlr0
        hdr['CDELT1'] = delt
        hdr['CTYPE1'] = 'WAVELENGTH'
        hdr['CRPIX2'] = 1
        hdr['CRVAL2'] = 1
        hdr['CDELT2'] = 1
        hdr['CTYPE2'] = 'PIXEL'
        return hdr

    def getHeaderList(self, image_list):
        final_list = []
        for image in image_list:
            for elem in image:
                if issubclass(fits.ImageHDU, type(elem)):
                    final_list.append(elem)
        return final_list

    @classmethod
    def types_getter(cls):
        from megaradrp.types import MasterBias, MasterDark, MasterBPM, MasterSlitFlat
        imgtypes = [None, MasterBPM, MasterBias, MasterDark, MasterSlitFlat]
        getters = [[cor.get_corrector_o, cor.get_corrector_t],
                   cor.get_corrector_p, cor.get_corrector_bias,
                   cor.get_corrector_dark, cor.get_corrector_sf
                   ]
        return imgtypes, getters

    @classmethod
    def load_getters(cls):
        import collections
        imgtypes, getters = cls.types_getter()
        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            print('load_getters', rtype, getter)
            if rtype is None:
                # Unconditional
                if isinstance(getter, collections.Iterable):
                    used_getters.extend(getter)
                else:
                    used_getters.append(getter)
            else:
                # Search
                for key, val in cls.RecipeInput.stored().items():
                    if isinstance(val.type, rtype):
                        if isinstance(getter, collections.Iterable):
                            used_getters.extend(getter)
                        else:
                            used_getters.append(getter)
                        break
                else:
                    pass
        return used_getters


    @classmethod
    def init_filters_generic(cls, rinput, getters, ins):
        from numina.flow import SerialFlow
        # with BPM, bias, dark, flat and sky
        #if emirdrp.ext.gtc.RUN_IN_GTC:
        #    _logger.debug('running in GTC environment')
        #else:
        cls.logger.debug('running outside of GTC environment')
        datamodel = MegaraDataModel()
        meta = datamodel.gather_info(rinput)
        cls.logger.debug('obresult info')
        for entry in meta['obresult']:
            cls.logger.debug('frame info is %s', entry)
        correctors = [getter(rinput, meta, ins) for getter in getters]

        flow = SerialFlow(correctors)

        return flow


    @classmethod
    def init_filters(cls, rinput, ins):
        getters = cls.load_getters()
        return cls.init_filters_generic(rinput, getters, ins)
