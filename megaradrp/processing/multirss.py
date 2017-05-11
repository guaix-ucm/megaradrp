#
# Copyright 2011-2017 Universidad Complutense de Madrid
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
