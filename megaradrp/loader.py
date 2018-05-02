#
# Copyright 2014-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Load MEGARA DRP"""

from numina.core import drp_load


def load_drp():
    """Entry point to load MEGARA DRP."""
    return drp_load('megaradrp', 'drp.yaml')
