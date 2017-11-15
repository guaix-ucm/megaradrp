#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""


import numina.types.structured as structured


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
