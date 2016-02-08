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

from ..core.processing import trim_and_o_hdu

_logger = logging.getLogger('megara.processing')


class OverscanCorrector(TagOptionalCorrector):

    '''A Node that corrects a frame from overscan.'''

    def __init__(self, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        # FIXME: these should come from the header
        bng = [1, 1]
        nr = 2056 / bng[0]
        nc = 2048 / bng[1]
        nr2 = 2 * nr
        nc2 = 2 * nc
        oscan1 = 50 / bng[0]
        oscan2 = oscan1 * 2
        psc1 = 50 / bng[0]
        psc2 = 2 * psc1
        fshape = (nr2 + oscan2, nc2 + psc2)
        # Row block 1
        rb1 = slice(0, nr)
        rb1m = slice(nr, nr + oscan1)
        # Row block 2
        rb2 = slice(nr + oscan2, nr2 + oscan2)
        rb2m = slice(nr + oscan1, nr + oscan2)
        # Col block
        cb = slice(psc1, nc2 + psc1)

        # Col block left
        cbl = slice(0, psc1)
        # Col block right
        cbr = slice(nc2 + psc1, nc2 + psc2)

        # Mode normal
        self.trim1 = (rb1, cb)
        self.pcol1 = (rb1, cbl)
        self.ocol1 = (rb1, cbr)
        self.orow1 = (rb1m, cb)

        self.trim2 = (rb2, cb)
        self.pcol2 = (rb2, cbr)
        self.ocol2 = (rb2, cbl)
        self.orow2 = (rb2m, cb)

        if tagger is None:
            tagger = TagFits('NUM-OVPE', 'Over scan/prescan')

        super(OverscanCorrector, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)

    def _run(self, img):
        data = img[0].data

        p1 = data[self.pcol1].mean()
        _logger.debug('prescan1 is %f', p1)
        or1 = data[self.orow1].mean()
        _logger.debug('row overscan1 is %f', or1)
        oc1 = data[self.ocol1].mean()
        _logger.debug('col overscan1 is %f', oc1)
        avg = (p1 + or1 + oc1) / 3.0
        _logger.debug('average scan1 is %f', avg)
        data[self.trim1] -= avg

        p2 = data[self.pcol2].mean()
        _logger.debug('prescan2 is %f', p2)
        or2 = data[self.orow2].mean()
        _logger.debug('row overscan2 is %f', or2)
        oc2 = data[self.ocol2].mean()
        _logger.debug('col overscan2 is %f', oc2)
        avg = (p2 + or2 + oc2) / 3.0
        _logger.debug('average scan2 is %f', avg)
        data[self.trim2] -= avg
        return img


class TrimImage(TagOptionalCorrector):

    '''A Node that trims images.'''

    def __init__(self, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-TRIM', 'Trimming')

        super(TrimImage, self).__init__(datamodel=datamodel,
                                        tagger=tagger,
                                        mark=mark,
                                        dtype=dtype)

    def _run(self, img):
        _logger.debug('trimming image %s', img)

        img[0] = trim_and_o_hdu(img[0])

        return img

