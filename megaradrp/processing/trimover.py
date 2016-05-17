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

import logging
import numpy as np
from numina.flow.processing import TagOptionalCorrector, TagFits

from megaradrp.core.processing import trimOut
from numina.core.config import  get_conf_value

_logger = logging.getLogger('megara.processing')


class OverscanCorrector(TagOptionalCorrector):
    '''A Node that corrects a frame from overscan.'''

    def __init__(self, datamodel=None, mark=True, tagger=None,
                 dtype='float32', confFile={}):

        trim1 = self._get_conf_value(confFile, 'trim1')
        trim2 = self._get_conf_value(confFile, 'trim2')
        bng = self._get_conf_value(confFile, 'bng')
        overscan1 = self._get_conf_value(confFile, 'overscan1')
        overscan2 = self._get_conf_value(confFile, 'overscan2')
        # psc1 = get_conf_value('MEGARA','default','prescan')
        prescan1 = self._get_conf_value(confFile, 'prescan1')
        prescan2 = self._get_conf_value(confFile, 'prescan2')
        middle1 = self._get_conf_value(confFile, 'middle1')
        middle2 = self._get_conf_value(confFile, 'middle2')

        bng = bng if bng else [1,1]
        trim1 = trim1 if trim1 else [[0,2056],[50,4146]]
        trim2 = trim2 if trim2 else [[2156,4212],[50,4146]]
        prescan1 = prescan1 if prescan1 else [[0,2056],[4146,4196]]
        prescan2 = prescan2 if prescan2 else [[2156,4212],[0,50]]
        overscan1 = overscan1 if overscan1 else [[0,2056],[0,50]]
        overscan2 = overscan2 if overscan2 else [[2156,4212],[4146,4196]]
        middle1 = middle1 if middle1 else [[2056,2106],[50,4146]]
        middle2 = middle2 if middle2 else [[2106,2156],[50,4146]]

        auxX,auxY,auxZ,auxT = self.data_binning(trim1,bng)
        middleX,middleY,middleZ,middleT = self.data_binning(middle1,bng)
        prescanX,prescanY,prescanZ,prescanT = self.data_binning(prescan1,bng)
        overscanX,overscanY,overscanZ,overscanT = self.data_binning(overscan1,bng)
        self.trim1 = (slice(auxX,auxY), slice(auxZ, auxT))
        self.orow1 = (slice(middleX,middleY), slice(middleZ, middleT))
        self.pcol1 = (slice(prescanX,prescanY), slice(prescanZ, prescanT))
        self.ocol1 = (slice(overscanX,overscanY), slice(overscanZ, overscanT))

        auxX,auxY,auxZ,auxT = self.data_binning(trim2,bng)
        middleX,middleY,middleZ,middleT = self.data_binning(middle2,bng)
        prescanX,prescanY,prescanZ,prescanT = self.data_binning(prescan2,bng)
        overscanX,overscanY,overscanZ,overscanT = self.data_binning(overscan2,bng)
        self.trim2 = (slice(auxX,auxY), slice(auxZ, auxT))
        self.orow2 = (slice(middleX,middleY), slice(middleZ, middleT))
        self.pcol2 = (slice(prescanX,prescanY), slice(prescanZ, prescanT))
        self.ocol2 = (slice(overscanX,overscanY), slice(overscanZ, overscanT))

        self.test_image()

        if tagger is None:
            tagger = TagFits('NUM-OVPE', 'Over scan/prescan')

        super(OverscanCorrector, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)

    def _get_conf_value(self, confFile, key=''):
        if confFile:
            if key in confFile.keys():
                return confFile[key]
            else:
                raise ValueError('Key is not in configuration file')
        raise ValueError('Instrument is not in the system')

    def data_binning(self, data, binning):
        '''
         Axis:x --> factorX
         Axis:y --> factorY
        '''
        factorX = 1.0/binning[1]
        factorY = 1.0/binning[0]
        x = int(factorY * data[0][0])
        y = int(factorX * data[0][1])
        z = int(factorY * data[1][0])
        t = int(factorX * data[1][1])

        return x,y,z,t

    def test_image(self):
        import astropy.io.fits as fits

        data = np.empty((4212,4196), dtype='float32')
        data[self.pcol1]+=1
        data[self.orow1]+=10
        data[self.ocol1]+=100
        data[self.trim1]+=1000

        data[self.pcol2]+=5
        data[self.orow2]+=50
        data[self.ocol2]+=500
        data[self.trim2]+=5000

        fits.writeto('eq_estimado.fits', data, clobber=True)

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

        super(TrimImage, self).__init__(datamodel=datamodel, tagger=tagger,
                                        mark=mark, dtype=dtype)

    def _run(self, img):
        _logger.debug('trimming image %s', img)

        img[0] = trimOut(img[0])

        return img
