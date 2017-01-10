#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

"""Products of the Megara Pipeline"""


from numina.core import DataFrameType, DataProductType
from numina.core.types import DataType
from numina.core.products import DataProductTag


class MEGARAProductFrame(DataProductTag, DataFrameType):
    """A MEGARA product image"""
    pass


class ProcessedFrame(DataFrameType):
    """A processed image not to be stored"""
    pass


class ProcessedRSS(ProcessedFrame):
    """A processed RSS image not to be stored"""
    pass


class ProcessedMultiRSS(ProcessedFrame):
    """A processed RSS image not to be stored"""
    pass


class MasterBias(MEGARAProductFrame):
    """A Master Bias image"""
    pass


class MasterTwilightFlat(MEGARAProductFrame):
    pass


class MasterDark(MEGARAProductFrame):
    """A Master Dark image"""
    pass


class MasterFiberFlat(MEGARAProductFrame):
    pass


class MasterSlitFlat(MEGARAProductFrame):
    pass


class MasterFiberFlatFrame(MEGARAProductFrame):
    pass


class MasterBPM(MEGARAProductFrame):
    """Bad Pixel Mask product"""
    pass


class MasterSensitivity(MEGARAProductFrame):
    pass


class MasterWeights(DataProductType):
    def __init__(self, default=None):
        super(MasterWeights, self).__init__(ptype=dict, default=default)

    def _datatype_dump(self, obj, where):
        import shutil

        filename = where.destination + '.tar'

        shutil.copy(obj, filename)
        shutil.rmtree(obj.split(filename)[0])
        return filename

    def _datatype_load(self, obj):
        try:
            import tarfile
            return tarfile.open(obj, 'r')
        except IOError as e:
            raise e


class JSONstorage(DataType):
    def __init__(self, default=None):
        super(JSONstorage, self).__init__(ptype=dict, default=default)

    def _datatype_dump(self, obj, where):
        import json
        filename = where.destination + '.json'

        with open(filename, 'w') as fd:
            fd.write(json.dumps(obj, sort_keys=True, indent=2,
                                separators=(',', ': ')))

        return filename

    def _datatype_load(self, obj):
        import json
        try:
            with open(obj, 'r') as fd:
                data = json.load(fd)
        except IOError as e:
            raise e
        return data


class LCBCalibration(JSONstorage):
    pass