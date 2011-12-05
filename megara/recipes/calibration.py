#
# Copyright 2011 Universidad Complutense de Madrid
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

'''Calibration Recipes for Megara'''

import logging

import numpy
import pyfits
from numina import RecipeBase, __version__
from numina.recipes import Parameter, log_to_history

from megara.products import MasterBias, MasterDark, MasterFlat

__all__ = ['BiasRecipe', 'DarkRecipe', 'FlatRecipe']

_logger = logging.getLogger('megara.recipes')

class BiasRecipe(RecipeBase):
    '''Process BIAS images and create MASTER_BIAS.'''

    __requires__ = []
    __provides__ = [MasterBias]

    def __init__(self):
        super(BiasRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    @log_to_history(_logger)
    def run(self, rb):
    	_logger.info('starting bias reduction')

        images = rb.images

        cdata = []

        try:
            for image in images:
                hdulist = pyfits.open(image, memmap=True, mode='readonly')
                cdata.append(hdulist)

            _logger.info('stacking images')
            data = numpy.zeros(cdata[0][1].data.shape, dtype='float32')
            for hdulist in cdata:
                data += hdulist[1].data

            data /= len(cdata)
            data += 2.0

            hdu = pyfits.ImageHDU(data, header=cdata[0][1].header)
    
            # update hdu header with
            # reduction keywords
            hdr = cdata[0][0].header
            hdr.update('FILENAME', 'master_bias-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'BIAS', 'Image type')
            hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'BiasRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([cdata[0][0], hdu])

            _logger.info('bias reduction ended')

            return {'products': [MasterBias(hdulist)]}
        finally:
            for hdulist in cdata:
                hdulist.close()

class DarkRecipe(RecipeBase):
    '''Process DARK images and provide MASTER_DARK. '''

    __requires__ = [Parameter('master_bias', MasterBias, 'comment')]
    __provides__ = [MasterDark]

    def __init__(self):
        super(DarkRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    @log_to_history(_logger)
    def run(self, block):

    	_logger.info('starting dark reduction')

        try:
            _logger.info('subtracting bias %s', str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd[1].data
                        data -= master_bias[1].data
                

            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd[1].data.copy()
                hdr = fd[1].header
                bhdr = fd[0].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd[1].data
                    data += add_data

            hdu = pyfits.ImageHDU(data, header=hdr)
    
            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_dark-%(block_id)d.fits' % self.environ)
            bhdr.update('FILENAME', 'master_dark-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'DARK', 'Image type')
            hdr.update('NUMTYP', 'MASTER_DARK', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'DarkRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=bhdr), hdu])

            _logger.info('dark reduction ended')



            return {'products': [MasterDark(hdulist)]}
        finally:
            pass

class FlatRecipe(RecipeBase):
    '''Process FLAT images and provide MASTER_FLAT. '''

    __requires__ = [
                    Parameter('master_bias', MasterBias, 'comment'),
                    Parameter('master_dark', MasterDark, 'comment')
                    ]
    __provides__ = [MasterFlat]

    def __init__(self):
        super(FlatRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    @log_to_history(_logger)
    def run(self, block):

    	_logger.info('starting flat reduction')

        try:
            _logger.info('subtracting bias %s', str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_bias['primary'].data
                

            _logger.info('subtracting dark %s', str(self.parameters['master_dark']))
            with pyfits.open(self.parameters['master_dark'], mode='readonly') as master_dark:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_dark['primary'].data


            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd['PRIMARY'].data.copy()
                hdr = fd['PRIMARY'].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd['primary'].data
                    data += add_data

            # Normalize flat to mean 1.0
            data[:] = 1.0

            hdu = pyfits.PrimaryHDU(data, header=hdr)

            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_flat-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'FLAT', 'Image type')
            hdr.update('NUMTYP', 'MASTER_FLAT', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'FlatRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('flat reduction ended')

            return {'products': [MasterFlat(hdulist)]}
        finally:
            pass

