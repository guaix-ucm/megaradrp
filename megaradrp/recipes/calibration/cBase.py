import logging

from astropy.io import fits


from numina.array.combine import median as c_median
from numina.core import BaseRecipeAutoQC

_logger = logging.getLogger('numina.recipes.megara')

class MegaraBaseRecipe (BaseRecipeAutoQC):

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