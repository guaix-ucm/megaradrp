#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy
import astropy.io.fits as fits

from megaradrp.utils import add_collapsed_mos_extension

def test_add_add_collapsed_mos_extension():

    img = fits.HDUList()

    data = numpy.ones((56, 1945), dtype='int16')
    for i in range(56):
        data[i] = i

    img.append(fits.ImageHDU(data))
    img.append(fits.ImageHDU(name='FIBERS'))

    result = add_collapsed_mos_extension(img)

    assert img[0].shape == (56, 1945)
    assert 'FIBERS' in result
    assert 'COLLAPSED' in result
    ext = result['COLLAPSED'].data
    assert ext.shape == (8, 1945)
    assert numpy.allclose(ext[0, :], 3.0)
    assert numpy.allclose(ext[1, :], 10.0)
    assert numpy.allclose(ext[2, :], 17.0)
    assert numpy.allclose(ext[3, :], 24.0)
    assert numpy.allclose(ext[4, :], 31.0)
    assert numpy.allclose(ext[5, :], 38.0)
    assert numpy.allclose(ext[6, :], 45.0)
    assert numpy.allclose(ext[7, :], 52.0)
