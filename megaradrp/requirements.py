#
# Copyright 2015-2016 Universidad Complutense de Madrid
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

"""Typical requirements of recipes"""

from numina.core import Requirement

import megaradrp.types
import megaradrp.products


class MasterBiasRequirement(Requirement):
    def __init__(self, optional=False):
        super(MasterBiasRequirement,
              self).__init__(megaradrp.types.MasterBias,
                             'Master BIAS image',
                             optional=optional
                             )


class MasterBPMRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBPMRequirement,
              self).__init__(megaradrp.types.MasterBPM,
                             'Master Bad Pixel Mask',
                             optional=optional
                             )


class MasterDarkRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterDarkRequirement,
              self).__init__(megaradrp.types.MasterDark, 'Master DARK image',
                             optional=optional)


class MasterFiberFlatRequirement(Requirement):
    def __init__(self):
        super(MasterFiberFlatRequirement,
              self).__init__(megaradrp.types.MasterFiberFlat,
                             'Master fiber flat calibration'
                             )


class MasterSlitFlatRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterSlitFlatRequirement,
              self).__init__(megaradrp.types.MasterSlitFlat,
                             'Master slit flat calibration',
                             optional=optional
                             )


class MasterFiberFlatFrameRequirement(Requirement):
    def __init__(self):
        super(MasterFiberFlatFrameRequirement,
              self).__init__(megaradrp.types.MasterFiberFlatFrame,'Master fiber flat frame')


class MasterTwilightRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterTwilightRequirement,
              self).__init__(megaradrp.types.MasterSlitFlat,'Master slit flat calibration')


class MasterTraceMapRequirement(Requirement):
    def __init__(self):
        super(MasterTraceMapRequirement,
              self).__init__(megaradrp.products.TraceMap, 'Trace information of the Apertures')
