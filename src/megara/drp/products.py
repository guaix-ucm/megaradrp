#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''Products of the Megara Pipeline'''


'''
    RAW_BIAS FrameDataProduct
    RAW_DARK FrameDataProduct
    RAW_FLAT FrameDataProduct
    RAW_ILLUM FrameDataProduct
    RAW_SCIENCE FrameDataProduct

    MASTER_BIAS  FrameDataProduct(detector)
    MASTER_DARK  FrameDataProduct(detector, exposure)
    MASTER_FLAT  FrameDataProduct(detector, grism)
    MASTER_ILLUM FrameDataProduct(detector, grism)

    POINTING FrameDataProduct
    MOSAIC FrameDataProduct

'''

import numpy

from numina.core import FrameDataProduct, DataProduct

class MasterBias(FrameDataProduct):
    pass

class MasterDark(FrameDataProduct):
    pass

class MasterFiberFlat(FrameDataProduct):
    pass

class MasterSensitivity(FrameDataProduct):
    pass

class TraceMapType(DataProduct):
    def __init__(self, default=None):
        super(TraceMapType, self).__init__(ptype=numpy.ndarray, default=default)

    def store(self, obj):
        return numpy.loadtxt(obj)
        


