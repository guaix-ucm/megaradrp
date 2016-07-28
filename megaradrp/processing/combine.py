#
# Copyright 2016 Universidad Complutense de Madrid
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


"""Combination routines"""

from __future__ import division
#
import logging

import numpy
from astropy.io import fits

from numina.array import combine


_logger = logging.getLogger('numina.recipes.megara')


def basic_processing(rinput, flow):

    cdata = []

    _logger.info('processing input images')
    for frame in rinput.obresult.images:
        hdulist = frame.open()
        fname = hdulist.filename()
        if fname:
            _logger.info('input is %s', fname)
        else:
            _logger.info('input is %s', hdulist)

        final = flow(hdulist)
        _logger.debug('output is input: %s', final is hdulist)

        cdata.append(final)

    return cdata


def basic_processing_with_combination(rinput, flow,
                                      method=combine.mean,
                                      errors=True):
    odata = []
    cdata = []
    try:
        _logger.info('processing input images')
        for frame in rinput.obresult.images:
            hdulist = frame.open()
            fname = hdulist.filename()
            if fname:
                _logger.info('input is %s', fname)
            else:
                _logger.info('input is %s', hdulist)

            final = flow(hdulist)
            _logger.debug('output is input: %s', final is hdulist)

            cdata.append(final)

            # Files to be closed at the end
            odata.append(hdulist)
            if final is not hdulist:
                odata.append(final)

        base_header = cdata[0][0].header.copy()
        _logger.info("stacking %d images using '%s'", len(cdata), method.func_name)
        data = method([d[0].data for d in cdata], dtype='float32')
        hdu = fits.PrimaryHDU(data[0], header=base_header)
        _logger.debug('update result header')
        hdu.header['history'] = "Combined %d images using '%s'" % (len(cdata), method.func_name)
        if errors:
            varhdu = fits.ImageHDU(data[1], name='VARIANCE')
            num = fits.ImageHDU(data[2], name='MAP')
            result = fits.HDUList([hdu, varhdu, num])
        else:
            result = fits.HDUList([hdu])
    finally:
        _logger.debug('closing images')
        for hdulist in odata:
            hdulist.close()

    return result
