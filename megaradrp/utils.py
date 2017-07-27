#
# Copyright 2017 Universidad Complutense de Madrid
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
        reshaped = img[0].data.reshape((-1, 7, 1, shape[1])).mean(axis=(1, 2))
    elif axis == 1:
        reshaped = img[0].data.reshape((shape[0], size, 1, -1)).mean(axis=(1, 2))
    else:
        raise ValueError('invalid axis={}'.format(axis))

    ext.data = reshaped

    img.append(ext)

    return img
