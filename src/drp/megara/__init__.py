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

'''The MEGARA Data Reduction Pipeline.'''

import logging
import pkgutil
import importlib
    
import yaml
from numina.core import BaseInstrument, BasePipeline, import_core

from megara import __version__

_modes = [i for i in yaml.load_all(pkgutil.get_data('megara',
                                                   'obsmodes.yaml'))
          if i.instrument == 'MEGARA']

_equiv_class = {}

for _m in _modes:
    _m.recipe_class = import_object(_m.recipe)
    _equiv_class[_m.key] = _m.recipe_class

del _m

class MegaraPipeline(BasePipeline):
    def __init__(self):
        super(MegaraPipeline, self).__init__(name='MEGARA', 
                version=__version__,
                recipes=_equiv_class)

class MEGARA_Instrument(BaseInstrument):
    name = 'MEGARA'
    modes = _modes
