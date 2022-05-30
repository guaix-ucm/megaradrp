#
# Copyright 2015-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import print_function

from astropy.io import fits
import numpy as np

# row / column
_binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
_direc = ['normal', 'mirror']


def create(image, direction='normal', bins='11'):
    """Create a image with overscan for testing."""

    if direction not in _direc:
        raise ValueError(f"{direction} must be either 'normal' or 'mirror'")

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = np.fliplr

    if bins not in _binning:
        raise ValueError(f"{bins} must be one if '11', '12', '21, '22'")

    bng = _binning[bins]

    nr = 2056 / bng[0]
    nc = 2048 / bng[1]

    nr2 = 2 * nr
    nc2 = 2 * nc

    oscan1 = 50 / bng[0]
    oscan2 = oscan1 * 2

    psc1 = 50 / bng[0]
    psc2 = 2 * psc1

    fshape = (nr2 + oscan2, nc2 + psc2)

    # Row block 1
    rb1 = slice(0, nr)
    rb1m = slice(nr, nr + oscan1)
    # Row block 2
    rb2 = slice(nr + oscan2, nr2 + oscan2)
    rb2m = slice(nr + oscan1, nr + oscan2)
    # Col block
    cb = slice(psc1, nc2 + psc1)
    # Col block left
    cbl = slice(0, psc1)
    # Col block right
    cbr = slice(nc2 + psc1, nc2 + psc2)

    # Mode normal
    trim1 = (rb1, cb)
    pcol1 = (rb1, cbl)
    ocol1 = (rb1, cbr)
    orow1 = (rb1m, cb)
    print(trim1, ocol1, orow1, pcol1)

    trim2 = (rb2, cb)
    pcol2 = (rb2, cbr)
    ocol2 = (rb2, cbl)
    orow2 = (rb2m, cb)
    print(trim2, ocol2, orow2, pcol2)

    finaldata = np.zeros(fshape, dtype='float32')

    finaldata[trim1] = direcfun(np.atleast_2d(np.arange(0, nc2)))
    finaldata[trim2] = direcfun(np.atleast_2d(np.arange(0, nc2)))

    finaldata[orow1] = 3
    finaldata[orow2] = 4

    finaldata[pcol1] = 5
    finaldata[pcol2] = 6

    finaldata[ocol1] = 7
    finaldata[ocol2] = 8

    hdu = fits.PrimaryHDU(data=finaldata)
    hdu.writeto(image, overwrite=True)
