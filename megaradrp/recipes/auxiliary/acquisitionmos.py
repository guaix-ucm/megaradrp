#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Acquire with MOS Recipe for Megara"""


from __future__ import division, print_function

import logging


from numina.core import RecipeError

from megaradrp.core.recipe import MegaraBaseRecipe


_logger = logging.getLogger('numina.recipes.megara')


class AcquireMOSRecipe(MegaraBaseRecipe):
    """Process Focus images and find best focus."""

    def __init__(self):
        super(AcquireMOSRecipe, self).__init__("0.1.0")

    def run(self, rinput):
        raise RecipeError('Recipe not implemented')
