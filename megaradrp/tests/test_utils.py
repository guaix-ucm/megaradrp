#
# Copyright 2015 Universidad Complutense de Madrid
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
