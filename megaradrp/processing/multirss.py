#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy
import astropy.io.fits as fits
import uuid


def generate_multi_rss(imgs):

    nimages = len(imgs)
    if nimages == 0:
        raise ValueError('need at least 1 image')
    refimg = imgs[0]

    dim0 = refimg[0].shape[0]
    multishape = (dim0 * nimages, refimg[0].shape[1])
    fdata = numpy.empty(multishape)

    for idx, img in enumerate(imgs):
        fdata[idx*dim0: (idx+1)*dim0] = img[0].data

    hdu = fits.PrimaryHDU(data=fdata, header=refimg[0].header)
    hdu.header['MEG-NRSS'] = nimages
    # Generate a new UUID
    hdu.header['UUID'] = str(uuid.uuid1())

    allhdus = [hdu]
    for idx, img in enumerate(imgs, 1):
        fibers = img['FIBERS']
        fibers.header['EXTNAME'] = 'FIBERS{}'.format(idx)
        allhdus.append(fibers)
    result = fits.HDUList(allhdus)
    return result
