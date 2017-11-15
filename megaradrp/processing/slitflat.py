#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging
import datetime

from .fiberflat import CommonFlatCorrector

_logger = logging.getLogger(__name__)


class SlitFlatCorrector(CommonFlatCorrector):
    """A Node that corrects a frame from slit flat."""

    def __init__(self, slitflat, datamodel=None, calibid='calibid-unknown',
                 dtype='float32'):

        super(SlitFlatCorrector, self).__init__(slitflat, datamodel, calibid, dtype)
        self.flattag = 'slitflat'

    def header_update(self, hdr, imgid):
        hdr['NUM-SLTF'] = self.calibid
        hdr['history'] = 'Slit flat correction {}'.format(imgid)
        hdr['history'] = 'Slit flat correction time {}'.format(datetime.datetime.utcnow().isoformat())
