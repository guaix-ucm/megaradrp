#
# Copyright 2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

import numpy


def compute_pa_from_ipa(ipa, ins_angle=-163.256):
    """Recompute the PA from IPA """
    # get IPA keyword
    pa = -ins_angle + ipa
    return pa


def update_wcs_from_ipa(hdr, pa):
    """Recompute the WCS rotations from PA"""
    pa_rad = numpy.deg2rad(pa)
    cos_pa = math.cos(pa_rad)
    sin_pa = math.sin(pa_rad)

    # Update PC_ keywords
    hdr['PC1_1'] = cos_pa
    hdr['PC2_2'] = cos_pa
    hdr['PC1_2'] = sin_pa
    hdr['PC2_1'] = -sin_pa
    # CDELT1 must be negative
    hdr['CDELT1'] = -abs(hdr['CDELT1'])

    return hdr