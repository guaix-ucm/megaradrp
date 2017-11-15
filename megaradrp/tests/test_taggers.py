#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numina.core

import astropy.io.fits as fits

from megaradrp.taggers import tagger_empty, tagger_vph
from megaradrp.taggers import tagger_base_image


def test_tagger_empty():

    img1 = fits.PrimaryHDU(data=[1, 2, 3], header=fits.Header())
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

    img1 = fits.PrimaryHDU(data=[1, 2, 3], header=fits.Header())
    img1.header['VPH'] = 'VPH405_LR'
    img1.header['READM'] = 'MOD1'
    frame1 = numina.core.DataFrame(frame=fits.HDUList(img1))

    img2 = fits.PrimaryHDU(data=[1, 2, 3], header=fits.Header())
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


def test_tagger_base_image():

    img1 = fits.PrimaryHDU(data=[1, 2, 3], header=fits.Header())
    img1.header['VPH'] = 'VPH405_LR'
    img1.header['INSMODE'] = 'MOS'
    img1.header['READM'] = 'MOD1'
    frame1 = numina.core.DataFrame(frame=fits.HDUList(img1))

    img2 = fits.PrimaryHDU(data=[1, 2, 3], header=fits.Header())
    img2.header['VPH'] = 'VPH405_LR'
    img2.header['INSMODE'] = 'MOS'
    img2.header['READM'] = 'MOD1'
    frame2 = numina.core.DataFrame(frame=fits.HDUList(img2))

    ob = numina.core.ObservationResult()
    ob.frames = [frame1, frame2]

    tags = tagger_base_image(ob)

    assert tags == {'vph': 'VPH405_LR', 'insmode': 'MOS'}


def test_tagger_base_image_no_images():

    ob = numina.core.ObservationResult()

    tags = tagger_base_image(ob)

    assert tags == {}

if __name__ == "__main__":
    test_tagger_vph()
