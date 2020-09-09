#
# Copyright 2017-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Fix RSS imperfections
"""
import math

import numpy as np

from megaradrp.instrument.focalplane import FocalPlaneConf
import megaradrp.processing.wcs as mwcs


def fix_missing_fiber(rss, fibid):
    """Interpolate missing fiber fibid"""
    hdr = rss['FIBERS'].header
    fp = FocalPlaneConf.from_img(rss)
    # Fibers around 623 are
    idxs = fp.nearby_fibers(fibid)
    avg = np.zeros_like(rss[0].data[fibid - 1])
    for idx in idxs:
        avg += rss[0].data[idx - 1]
    avg /= len(idxs)

    l1fmt = "FIB{:03d}W1"
    l2fmt = "FIB{:03d}W2"
    try:
        # Change limits in header to the min-max of surrounding fibers
        for sub, func in [(l1fmt, max), (l2fmt, min)]:
            keys = []
            for idx in idxs:
                keys.append(hdr[sub.format(idx)])
            hdr[sub.format(fibid)] = func(keys)

        l1 = hdr[l1fmt.format(fibid)]
        l2 = hdr[l2fmt.format(fibid)]
        if l1 > 1:
            avg[0:l1 - 1] = 0
        if l2 < len(avg):
            avg[l2 - 1:] = 0
    except KeyError:
        pass

    rss[0].data[fibid - 1] = avg
    return rss


def recompute_wcs(hdr, ipa):
    """Recompute the WCS rotations from IPA """
    pa = mwcs.compute_pa_from_ipa(ipa)
    # print('IPA angle is:', ipa, 'PA angle is', math.fmod(pa, 360))
    x = hdr['PC1_1']
    y = hdr['PC1_2']
    # print('PA from header is:', np.rad2deg(math.atan2(y, x)))
    return mwcs.update_wcs_from_ipa(hdr, pa)
