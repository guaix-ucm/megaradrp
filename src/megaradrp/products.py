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
    RAW_BIAS DataFrameType
    RAW_DARK DataFrameType
    RAW_FLAT DataFrameType
    RAW_ILLUM DataFrameType
    RAW_SCIENCE DataFrameType

    MASTER_BIAS  DataFrameType(detector)
    MASTER_DARK  DataFrameType(detector, exposure)
    MASTER_FLAT  DataFrameType(detector, grism)
    MASTER_ILLUM DataFrameType(detector, grism)

    POINTING DataFrameType
    MOSAIC DataFrameType

'''


from numina.core import DataFrameType, DataProductType


class MasterBias(DataFrameType):
    pass


class MasterDark(DataFrameType):
    pass


class MasterFiberFlat(DataFrameType):
    pass


class MasterSensitivity(DataFrameType):
    pass


class TraceMapType(DataProductType):

    def __init__(self, default=None):
        super(TraceMapType, self).__init__(
            ptype=dict, default=default)
