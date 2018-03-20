#
# Copyright 2015-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import pytest

import numpy
import astropy.io.fits as fits

from megaradrp.utils import add_collapsed_mos_extension


def test_add_add_collapsed_mos_extension0():

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
    for i in range(7):
        assert numpy.allclose(ext[i, :], 3.0 + 7 * i)


def test_add_add_collapsed_mos_extension1():

    img = fits.HDUList()
    size0 = 1945

    size1 = 56
    size = 7
    ncol = size1 // size

    data = numpy.ones((size0, size1), dtype='int16')
    for i in range(size1):
        data[:, i] = i

    img.append(fits.ImageHDU(data))
    img.append(fits.ImageHDU(name='FIBERS'))

    result = add_collapsed_mos_extension(img, axis=1)

    assert img[0].shape == (size0, size1)
    assert 'FIBERS' in result
    assert 'COLLAPSED' in result
    ext = result['COLLAPSED'].data
    assert ext.shape == (size0, ncol)
    for i in range(size):
        assert numpy.allclose(ext[:, i], 24.0 + i)


@pytest.mark.parametrize("axis", [-1, 2])
def test_add_add_collapsed_mos_extension2(axis):

    img = fits.HDUList()

    data = numpy.ones((56, 1945), dtype='int16')
    for i in range(56):
        data[i] = i

    img.append(fits.ImageHDU(data))
    img.append(fits.ImageHDU(name='FIBERS'))

    with pytest.raises(ValueError):
        add_collapsed_mos_extension(img, axis=axis)
