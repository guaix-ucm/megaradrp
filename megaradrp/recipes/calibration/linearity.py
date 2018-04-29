#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Calibrate the non linearity correciont of MEGARA Detector"""


from megaradrp.core.recipe import MegaraBaseRecipe


class LinearityTestRecipe(MegaraBaseRecipe):
    """Process LINTEST images and create LINEARITY_CORRECTON."""
    pass


