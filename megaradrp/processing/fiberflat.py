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
import numpy as np

from numina.flow.processing import TagOptionalCorrector, TagFits


_logger = logging.getLogger('numina.processing')


class FiberFlatCorrector(TagOptionalCorrector):

    '''A Node that corrects from fiber flat.'''

    def __init__(self, fiberflat, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-MFF', 'MEGARA Fiber flat correction')

        super(FiberFlatCorrector, self).__init__(datamodel=datamodel,
                                                 tagger=tagger,
                                                 mark=mark,
                                                 dtype=dtype)

        if isinstance(fiberflat, fits.HDUList):
            self.corr = fiberflat[0].data
        elif isinstance(fiberflat, np.ndarray):
            self.corr = fiberflat
        self.corrmean = self.corr.mean()
        self.corrid = self.get_imgid(fiberflat)

    def _run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct from fiber flat in image %s', imgid)

        img[0].data = img[0].data / self.corr

        return img

