#
# Copyright 2016-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Combination routines"""

from __future__ import division

import datetime
import logging
import uuid

from astropy.io import fits

from numina.array import combine

import megaradrp.datamodel


def basic_processing_with_combination(rinput, flow,
                                      method=combine.mean,
                                      errors=True,
                                      prolog=None):
    return basic_processing_with_combination_frames(rinput.obresult.frames,
                                                    flow, method=method,
                                                    errors=errors,
                                                    prolog=prolog)


def basic_processing_with_combination_frames(frames,
                                             flow,
                                             method=combine.mean,
                                             errors=True,
                                             prolog=None
                                             ):

    _logger = logging.getLogger(__name__)
    odata = []
    cdata = []

    datamodel = megaradrp.datamodel.MegaraDataModel()
    try:
        _logger.info('processing input images')
        for frame in frames:
            hdulist = frame.open()
            fname = datamodel.get_imgid(hdulist)
            _logger.info('input is %s', fname)
            final = flow(hdulist)
            _logger.debug('output is input: %s', final is hdulist)
            cdata.append(final)
            # Files to be closed at the end
            odata.append(hdulist)
            if final is not hdulist:
                odata.append(final)

        first_image = cdata[0]
        base_header = first_image[0].header.copy()
        cnum = len(cdata)
        _logger.info("stacking %d images using '%s'", cnum, method.__name__)
        data = method([d[0].data for d in cdata], dtype='float32')
        hdu = fits.PrimaryHDU(data[0], header=base_header)
        _logger.debug('update result header')
        if prolog:
            _logger.debug('write prolog')
            hdu.header['history'] = prolog
        hdu.header['history'] = "Combined %d images using '%s'" % (cnum, method.__name__)
        hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())
        for img in cdata:
            hdu.header['history'] = "Image {}".format(datamodel.get_imgid(img))
        prevnum = base_header.get('NUM-NCOM', 1)
        hdu.header['NUM-NCOM'] = prevnum * cnum
        hdu.header['UUID'] = str(uuid.uuid1())

        # Copy extensions and then append 'variance' and 'map'
        result = fits.HDUList([hdu])
        for hdu in first_image[1:]:
            result.append(hdu.copy())

        # Headers of last image
#        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        # Append error extensions
        if errors:
            varhdu = fits.ImageHDU(data[1], name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(data[2].astype('int16'), name='MAP')
            result.append(num)

    finally:
        _logger.debug('closing images')
        for hdulist in odata:
            hdulist.close()

    return result
