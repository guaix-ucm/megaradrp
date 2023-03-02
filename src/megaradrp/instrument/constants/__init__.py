#
# Copyright 2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""MEGARA constants with units"""

import astropy.units as u


# FIXME: duplicated in megaradrp.instrument
# without units

# Platescale in focal plane of Folded-Cass
GTC_FC_A_PLATESCALE = 1.212 * u.arcsec / u.mm

# Reference instrument aligment angle
MEGARA_IAA = -163.854 * u.deg

# mm from center to center, upwards
SPAXEL_SCALE = 0.443  * u.mm