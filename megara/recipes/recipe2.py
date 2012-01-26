#
# Copyright 2011-2012 Universidad Complutense de Madrid
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

from numina import RecipeBase

__all__ = ['Recipe']

class Recipe(RecipeBase):
    def __init__(self, pp, cp):
        pass

    @classmethod
    def requires(cls):
        return []

    def run(self, rb):
        return {'result': {'direct_image': 0, 'qa': 1}}

