#
# Copyright 2018-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import math

import numpy

import megaradrp.instrument.constants as cons


def compute_pa_from_ipa(ipa, iaa=cons.MEGARA_IAA.value):
    """Recompute the PA from IPA

    Parameters
    ==========
    ipa: float
        Instrument Position Angle
    iaa: float
        Instrument Alignment Angle

    Returns
    =======
    position angle

    Compute the Position angle of the image from the IPA and
    the (generally fixed) Instrument Alignment Angle.
    The angles must be in the same units.

    """

    pa = -iaa + ipa
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