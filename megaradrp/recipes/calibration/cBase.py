import logging

from astropy.io import fits

from megaradrp.processing import OverscanCorrector, TrimImage

from numina.array.combine import median as c_median
from numina.core import BaseRecipeAutoQC
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

_logger = logging.getLogger('numina.recipes.megara')

class MegaraBaseRecipe (BaseRecipeAutoQC):

    def __init__(self, version):
        self.__flow = {'BiasRecipe':[OverscanCorrector, TrimImage],
                       'ArcCalibrationRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       'FiberFlatRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       'TraceMapRecipe':[OverscanCorrector, TrimImage, BiasCorrector],
                       }
        super(MegaraBaseRecipe, self).__init__(version=version)

    def __generate_flow(self, params):
        flow = self.__flow[self.__class__.__name__]
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

        return basicflow

    def bias_process_common(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')

        with master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        hdu, data = self.hdu_creation(obresult, {'biasmap':mbias})

        hdr = hdu.header
        hdr['IMGTYP'] = ('FIBER_FLAT', 'Image type')
        hdr['NUMTYP'] = ('MASTER_FIBER_FLAT', 'Data product type')
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        reduced = fits.HDUList([hdu, varhdu, num])

        _logger.info('prereduction ended')
        _logger.info('extract fibers')

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

            data = c_median([d[0].data for d in cdata], dtype='float32')
            template_header = cdata[0][0].header
            hdu = fits.PrimaryHDU(data[0], header=template_header)
        finally:
            for hdulist in cdata:
                hdulist.close()

        return hdu, data