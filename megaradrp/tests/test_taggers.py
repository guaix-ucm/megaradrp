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

import numina.core

import astropy.io.fits as fits

from megaradrp.taggers import tagger_empty, tagger_vph

def test_tagger_empty():

    img1 = fits.PrimaryHDU(data=[1,2,3], header=fits.Header())
    img1.header['FILTER'] = 'FILTER-A'
    img1.header['READM'] = 'MOD1'
    frame1 = numina.core.DataFrame(frame=fits.HDUList(img1))

    ob = numina.core.ObservationResult()
    ob.frames = [frame1]

    tags = tagger_empty(ob)

    assert tags == {}


def test_tagger_empty_no_images():

    ob = numina.core.ObservationResult()

    tags = tagger_empty(ob)

    assert tags == {}


def test_tagger_vph():

    img1 = fits.PrimaryHDU(data=[1,2,3], header=fits.Header())
    img1.header['VPH'] = 'VPH405_LR'
    img1.header['READM'] = 'MOD1'
    frame1 = numina.core.DataFrame(frame=fits.HDUList(img1))

    img2 = fits.PrimaryHDU(data=[1,2,3], header=fits.Header())
    img2.header['VPH'] = 'VPH405_LR'
    img2.header['READM'] = 'MOD1'
    frame2 = numina.core.DataFrame(frame=fits.HDUList(img2))

    ob = numina.core.ObservationResult()
    ob.frames = [frame1, frame2]

    tags = tagger_vph(ob)

    assert tags == {'vph': 'VPH405_LR'}


def test_tagger_vph_no_images():

    ob = numina.core.ObservationResult()

    tags = tagger_vph(ob)

    assert tags == {}


if __name__ == "__main__":
    test_tagger_vph()

