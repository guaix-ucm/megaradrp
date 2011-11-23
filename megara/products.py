#
# Copyright 2011 Sergio Pascual
# 
# This file is part of Pontifex
# 
# Pontifex is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Pontifex is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Pontifex.  If not, see <http://www.gnu.org/licenses/>.
#

'''Products of the Megara Pipeline'''


'''
    RAW_BIAS Image
    RAW_DARK Image
    RAW_FLAT Image
    RAW_ILLUM Image
    RAW_SCIENCE Image

    MASTER_BIAS  Image(detector)
    MASTER_DARK  Image(detector, exposure)
    MASTER_FLAT  Image(detector, grism)
    MASTER_ILLUM Image(detector, grism)

    POINTING Image
    MOSAIC Image

'''

from numina.recipes import Image

class MasterBias(Image):
    def __init__(self, hdu):
        super(MasterBias, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['readmode']

class MasterDark(Image):
    def __init__(self, hdu):
        super(MasterDark, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['readmode']

class MasterFlat(Image):
    def __init__(self, hdu):
        super(MasterFlat, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['readmode']
        yield 'grism0', hdr['grism']

class MasterIllum(Image):
    def __init__(self, hdu):
        super(MasterIllum, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['readmode']
        yield 'grism0', hdr['grism']

