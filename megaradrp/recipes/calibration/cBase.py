import logging

from astropy.io import fits

from numina.array.combine import median as c_median
from numina.core import BaseRecipe
from numina.core.dataholders import Product
from numina.core.metarecipes import RecipeType
from numina.core.products import QualityControlProduct
from numina.core.recipeinout import add_product

from six import with_metaclass

_logger = logging.getLogger('numina.recipes.megara')

@add_product(qc=Product(QualityControlProduct, dest='qc'))
class MegaraBaseRecipe (with_metaclass(RecipeType, BaseRecipe)):

    def __init__(self, version):
        super(MegaraBaseRecipe, self).__init__(version=version)

    def hdu_creation(self, obresult, basicflow):
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