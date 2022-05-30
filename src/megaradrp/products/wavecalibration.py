#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline: Wavelength  Calibration"""


from numina.array.wavecalib.arccalibration import SolutionArcCalibration
import numpy.polynomial.polynomial as nppol

from .structured import BaseStructuredCalibration


class FiberSolutionArcCalibration(object):
    def __init__(self, fibid, solution):
        self.fibid = fibid
        self.solution = solution

    def __getstate__(self):
        return {'fibid': self.fibid,
                'solution': self.solution.__getstate__()
        }

    def __setstate__(self, state):
        self.fibid = state['fibid']
        new = SolutionArcCalibration.__new__(SolutionArcCalibration)
        new.__setstate__(state['solution'])
        self.solution = new


class WavelengthCalibration(BaseStructuredCalibration):
    """Wavelength Calibration Product
    """

    __tags__ = ['insmode', 'vph']

    def __init__(self, instrument='unknown'):
        super(WavelengthCalibration, self).__init__(instrument)
        self.contents = []
        self.global_offset = nppol.Polynomial([0.0])

    def tag_names(self):
        return ['insmode', 'vph']

    def __getstate__(self):
        st = super(WavelengthCalibration, self).__getstate__()

        st['contents'] = [val.__getstate__() for val in self.contents]
        st['global_offset'] = self.global_offset.coef
        return st

    def __setstate__(self, state):
        super(WavelengthCalibration, self).__setstate__(state)

        self.global_offset = nppol.Polynomial(state.get('global_offset', [0.0]))

        self.contents = []
        # Handle dictionary
        if isinstance(state['contents'], dict):
            values = state['contents'].values()
        else:
            values = state['contents']

        for val in values:
            new = FiberSolutionArcCalibration.__new__(FiberSolutionArcCalibration)
            new.__setstate__(val)
            self.contents.append(new)
