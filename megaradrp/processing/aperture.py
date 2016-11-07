#
# Copyright 2011-2016 Universidad Complutense de Madrid
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
import datetime

from numina.flow.processing import Corrector

from megaradrp.core.processing import apextract_tracemap_2

_logger = logging.getLogger(__name__)


class ApertureExtractor(Corrector):
    """A Node that extracts apertures."""

    def __init__(self, tracemap, datamodel=None, dtype='float32'):

        self.tracemap = tracemap

        super(ApertureExtractor, self).__init__(datamodel=datamodel,
                                                calibid=tracemap.uuid,
                                                dtype=dtype)

    def run(self, img):
        _logger.debug('simple aperture extraction')
        imgid = self.get_imgid(img)
        _logger.debug('extracting (apextract_tracemap) in image %s', imgid)
        _logger.debug('with trace map %s', self.calibid)
        rssdata = apextract_tracemap_2(img[0].data, self.tracemap)
        img[0].data = rssdata

        hdr = img[0].header

        hdr['NUM-APE'] = self.calibid
        hdr['history'] = 'Aperture extraction with {}'.format(self.calibid)
        hdr['history'] = 'Aperture extraction time {}'.format(datetime.datetime.utcnow().isoformat())

        _logger.debug('remove VARIANCE and MAP extensions')
        del img['variance']
        del img['map']

        # Update Fibers
        fibers_ext = img['FIBERS']
        fibers_ext_headers = fibers_ext.header
        for aper in self.tracemap.contents:
            key = "FIB%03d_V" % aper.fibid
            fibers_ext_headers[key] = aper.valid
            key = "FIB%03dS1" % aper.fibid
            fibers_ext_headers[key] = aper.start
            key = "FIB%03dS2" % aper.fibid
            fibers_ext_headers[key] = aper.stop

        return img
