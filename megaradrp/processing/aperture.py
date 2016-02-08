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

from numina.flow.processing import TagOptionalCorrector, TagFits

from megaradrp.core.processing import apextract, apextract_tracemap


_logger = logging.getLogger('numina.processing')


class ApertureExtractor(TagOptionalCorrector):
    """A Node that extracts apertures."""

    def __init__(self, trace, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-MAE', 'MEGARA Aperture extractor')

        super(ApertureExtractor, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)
        self.trace = trace

    def _run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('extracting apertures in image %s', imgid)
        rss = apextract(img[0].data, self.trace)
        img[0].data = rss

        return img


class ApertureExtractor2(TagOptionalCorrector):
    """A Node that extracts apertures."""

    def __init__(self, trace, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-MAE', 'MEGARA Aperture extractor')
        
        self.trace = trace

        super(ApertureExtractor2, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)


    def _run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('extracting (apextract_tracemap) in image %s', imgid)
        rss = apextract_tracemap(img[0].data, self.trace)
        img[0].data = rss
        return img

