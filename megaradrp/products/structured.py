#
# Copyright 2011-2017 Universidad Complutense de Madrid
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


import numina.core.products.structured as structured


class BaseStructuredCalibration(structured.BaseStructuredCalibration):
    def __init__(self, instrument='unknown'):
        super(BaseStructuredCalibration, self).__init__(instrument)
        self.total_fibers = 0
        self.missing_fibers = []
        self.error_fitting = []

    def __getstate__(self):
        st = super(BaseStructuredCalibration, self).__getstate__()

        keys = ["total_fibers", "missing_fibers", "error_fitting"]

        for key in keys:
            st[key] = self.__dict__[key]

        return st
