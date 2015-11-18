import logging

from astropy.io import fits

from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe
from megaradrp.products import MasterBias

from numina.core import Product, RecipeError
from numina.core.requirements import ObservationResultRequirement

_logger = logging.getLogger('numina.recipes.megara')


class BiasRecipe(MegaraBaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    obresult = ObservationResultRequirement()
    biasframe = Product(MasterBias)

    def __init__(self):
        super(BiasRecipe, self).__init__(version="0.1.0")

    def run(self, rinput):

        _logger.info('starting bias reduction')

        if not rinput.obresult.images:
            raise RecipeError('Frame list is empty')

        hdu, data = self.hdu_creation(rinput.obresult)

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