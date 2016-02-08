#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

from astropy.io import fits
import numina.array.combine as combine
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector, BadPixelCorrector, \
    DarkCorrector
from numina.core import BaseRecipe
from numina.core.dataholders import Product
from numina.core.products import QualityControlProduct

from megaradrp.processing.trimover import OverscanCorrector, TrimImage

_logger = logging.getLogger('numina.recipes.megara')


class MegaraBaseRecipe(BaseRecipe):
    """Base clase for all MEGARA Recipes"""

    qc = Product(QualityControlProduct, dest='qc')

    def __init__(self, version):
        self.__flow = {'ArcCalibrationRecipe': [OverscanCorrector, TrimImage,
                                                BiasCorrector,
                                                BadPixelCorrector,
                                                DarkCorrector],
                       'BadPixelsMaskRecipe': [OverscanCorrector, TrimImage,
                                               BiasCorrector, DarkCorrector],
                       'BiasRecipe': [OverscanCorrector, TrimImage,
                                      BadPixelCorrector],
                       'DarkRecipe': [OverscanCorrector, TrimImage,
                                      BiasCorrector],
                       'FiberFlatRecipe': [OverscanCorrector, TrimImage,
                                           BiasCorrector, BadPixelCorrector,
                                           DarkCorrector],
                       'TraceMapRecipe': [OverscanCorrector, TrimImage,
                                          BiasCorrector, BadPixelCorrector,
                                          DarkCorrector],

                       }
        super(MegaraBaseRecipe, self).__init__(version=version)

    def __generate_flow(self, params):
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
                elif issubclass(TrimImage, flow[cont]) or issubclass(
                        OverscanCorrector, flow[cont]):
                    flow[cont] = (flow[cont]())
                cont += 1
            basicflow = SerialFlow(flow)

        except Exception as e:
            _logger.error(e)
            raise (e)
        del flow
        return basicflow

    def bias_process_common(self, obresult, img):

        # with master_bias.open() as hdul:
        #     mbias = hdul[0].data.copy()

        hdu, data = self.hdu_creation(obresult, img)

        hdr = hdu[0].header
        # FIXME: this is incorrect in general
        hdr['IMGTYP'] = ('FIBER_FLAT', 'Image type')
        hdr['NUMTYP'] = ('MASTER_FIBER_FLAT', 'Data product type')
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        reduced = fits.HDUList(hdu + [varhdu, num])
        return reduced

    def hdu_creation(self, obresult, params={}):

        basicflow = self.__generate_flow(params)
        lista = []
        cdata = []
        try:
            for frame in obresult.images:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                cdata.append(hdulist)

            _logger.info('stacking %d images using median', len(cdata))

            data = combine.median([d[0].data for d in cdata], dtype='float32')
            template_header = cdata[0][0].header
            hdu = fits.PrimaryHDU(data[0], header=template_header)
            lista.append(hdu)
            # if has_mask:
            #     hdumask = fits.ImageHDU(themask)
            #     lista.append(hdumask)
        finally:
            for hdulist in cdata:
                hdulist.close()

        return fits.HDUList(lista), data

    def get_parameters(self, rinput):

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        parameters = {'biasmap':mbias}

        if rinput.master_bpm:
            with rinput.master_bpm.open() as hdul:
                parameters['bpm'] = hdul[0].data.copy()
        if rinput.master_dark:
            with rinput.master_dark.open() as hdul:
                parameters['dark'] = hdul[0].data.copy()

        return parameters
