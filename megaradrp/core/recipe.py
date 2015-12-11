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
from numina.flow.processing import BiasCorrector
from numina.core import BaseRecipe
from numina.core.dataholders import Product
from numina.core.products import QualityControlProduct

from megaradrp.processing import OverscanCorrector, TrimImage

_logger = logging.getLogger('numina.recipes.megara')

class MegaraBaseRecipe(BaseRecipe):
    """Base clase for all MEGARA Recipes"""

    qc = Product(QualityControlProduct, dest='qc')

    def __init__(self, version):
        self.__flow = {'BiasRecipe':[OverscanCorrector, TrimImage],
                       'ArcCalibrationRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       'FiberFlatRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       'TraceMapRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       'BadPixelsMaskRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       }
        super(MegaraBaseRecipe, self).__init__(version=version)

    def __generate_flow(self, params):
        import copy
        ff  = self.__flow[self.__class__.__name__]
        flow = copy.deepcopy(ff)
        # flow = self.__flow[self.__class__.__name__]
        try:
            for cont in range(len(flow)):
                if issubclass(BiasCorrector, flow[cont]):
                    flow[cont] = (flow[cont](params['biasmap']))
                elif issubclass(TrimImage, flow[cont]) or issubclass(OverscanCorrector, flow[cont]):
                    flow[cont] = (flow[cont]())
            basicflow = SerialFlow(flow)

        except Exception as e:
            _logger.error(e)
            raise(e)
        del flow
        return basicflow

    def bias_process_common(self, obresult, master_bias):

        with master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        hdu, data = self.hdu_creation(obresult, {'biasmap':mbias})

        hdr = hdu.header
        # FIXME: this is incorrect in general
        hdr['IMGTYP'] = ('FIBER_FLAT', 'Image type')
        hdr['NUMTYP'] = ('MASTER_FIBER_FLAT', 'Data product type')
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        reduced = fits.HDUList([hdu, varhdu, num])
        return reduced

    def hdu_creation(self, obresult, params=None):

        basicflow = self.__generate_flow(params)

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
        finally:
            for hdulist in cdata:
                hdulist.close()

        return hdu, data