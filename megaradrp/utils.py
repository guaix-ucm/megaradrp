#
# Copyright 2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Some utils"""


import astropy.io.fits as fits


def add_collapsed_mos_extension(img, size=7, axis=0):
    """Add a collapsed image extension

    Parameters
    ==========
    img: astropy.io.fits.HDUList

    Returns
    =======
    astropy.io.fits.HDUList
        Updated image

    """

    ext = fits.ImageHDU(header=img[0].header, name='COLLAPSED')

    shape = img[0].data.shape
    if axis == 0:
        reshaped = img[0].data.reshape((-1, size, 1, shape[1])).mean(axis=(1, 2))
    elif axis == 1:
        reshaped = img[0].data.reshape((shape[0], size, 1, -1)).mean(axis=(1, 2))
    else:
        raise ValueError('invalid axis={}'.format(axis))

    ext.data = reshaped

    img.append(ext)

    return img
