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

'''Products of the Megara Pipeline'''


'''
    RAW_BIAS DataFrame
    RAW_DARK DataFrame
    RAW_FLAT DataFrame
    RAW_ILLUM DataFrame
    RAW_SCIENCE DataFrame

    MASTER_BIAS  DataFrame(detector)
    MASTER_DARK  DataFrame(detector, exposure)
    MASTER_FLAT  DataFrame(detector, grism)
    MASTER_ILLUM DataFrame(detector, grism)

    POINTING DataFrame
    MOSAIC DataFrame

'''

from numina.recipes import DataFrame

class MasterBias(DataFrame):
    def __init__(self, hdu):
        super(MasterBias, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[1].header
        yield 'spec1.detector.mode', hdr['readmode']

class MasterDark(DataFrame):
    def __init__(self, hdu):
        super(MasterDark, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[1].header
        yield 'spec1.detector.mode', hdr['readmode']

class MasterFlat(DataFrame):
    def __init__(self, hdu):
        super(MasterFlat, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'spec1.detector.mode', hdr['readmode']
        yield 'spec1.grism', hdr['grism']

class MasterIllum(DataFrame):
    def __init__(self, hdu):
        super(MasterIllum, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[1].header
        yield 'spec1.detector.mode', hdr['readmode']
        yield 'spec1.grism', hdr['grism']

