__author__ = 'Pica4x6'

import logging

from astropy.io import fits

from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
from numina.array.combine import median as c_median
from numina.flow import SerialFlow

from cBase import BaseRecipeAutoQC as MegaraBaseRecipe
from megaradrp.core import OverscanCorrector, TrimImage
from megaradrp.products import MasterBias

_logger = logging.getLogger('numina.recipes.megara')

class BiasRecipe(MegaraBaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    obresult = ObservationResultRequirement()

    biasframe = Product(MasterBias)

    def __init__(self):
        super(BiasRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):

        obresult = rinput.obresult
        _logger.info('starting bias reduction')

        if not obresult.images:
            raise RecipeError('Frame list is empty')

        cdata = []

        o_c = OverscanCorrector()
        t_i = TrimImage()

        basicflow = SerialFlow([o_c, t_i])

        try:
            for frame in obresult.images:
                print 'FRame', frame
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

        hdr = hdu.header
        hdr = self.set_base_headers(hdr)
        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        hdulist = fits.HDUList([hdu, varhdu, num])
        _logger.info('bias reduction ended')

        result = self.create_result(biasframe=hdu)
        return result