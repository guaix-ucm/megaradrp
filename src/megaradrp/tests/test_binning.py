#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy as np

from megaradrp.instrument.components.detector import binning

def test_binning():
    nr = 6
    nc = 8
    bc = 2
    br = 2
    bnr = nr / br
    bnc = nc / bc
    arr = np.arange(nr*nc, dtype='int64').reshape(nr, nc)

    barr = binning(arr, bc, br)

    assert barr.shape == (bnr, bnc, br, bc)

    assert np.all(barr[0,0] == np.array([[0, 1], [8,9]]))

    carr = barr.reshape(barr.shape[0], barr.shape[1], -1)

    res = np.array([[ 18,  26,  34,  42], [ 82,  90,  98, 106], [146, 154, 162, 170]])

    assert np.all(carr.sum(axis=-1) == res)


if __name__ == "__main__":
    test_binning()