#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

from astropy.io import fits
import numpy as np

from numina.flow.processing import Corrector


_logger = logging.getLogger('numina.processing')


class TwilightCorrector(Corrector):

    '''A Node that corrects from twilight.'''

    def __init__(self, twilight, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        super(TwilightCorrector, self).__init__(datamodel=datamodel,
                                                 dtype=dtype)

        if isinstance(twilight, fits.HDUList):
            self.corr = twilight[0].data
        elif isinstance(twilight, np.ndarray):
            self.corr = twilight
        self.corrmean = self.corr.mean()
        # self.corrid = self.get_imgid(fiberflat)

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct from twilight in image %s', imgid)

        # Avoid nan values when divide
        my_mask = self.corr == 0.0
        self.corr[my_mask] = 1.0

        img[0].data /= self.corr

        return img
