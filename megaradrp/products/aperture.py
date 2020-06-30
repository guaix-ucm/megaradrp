#
# Copyright 2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""


class GeometricAperture(object):
    def __init__(self, fibid, boxid, start, stop):
        self.fibid = fibid
        self.boxid = boxid
        self.start = start
        self.stop = stop

    @property
    def valid(self):
        return self.is_valid()

    def aper_center(self):
        raise NotImplementedError

    def is_valid(self):
        return True

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__ = state
