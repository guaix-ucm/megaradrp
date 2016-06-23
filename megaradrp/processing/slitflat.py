import logging
from numina.flow.processing import Corrector
import numpy as np

_logger = logging.getLogger('numina.processing')

class SlitFlatCorrector(Corrector):
    '''A Node that corrects a frame from slit flat.'''

    def __init__(self, slitflat, datamodel=None,
                 dtype='float32'):

        super(SlitFlatCorrector, self).__init__(datamodel, dtype)

        self.slitflat = slitflat

    def run(self, img):
        imgid = self.get_imgid(img)

        _logger.debug('correcting slit flat in %s', imgid)

        # Avoid nan values when divide
        my_mask = self.slitflat == 0.0
        self.slitflat[my_mask] = 1.0

        img[0].data /= self.slitflat

        return img