#
# Copyright 2017-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""

import numpy.polynomial.polynomial as nppol

from numina.util.convertfunc import json_serial_function, convert_function

from megaradrp.processing.modelmap import calc_matrix_cols, aper_extract
from .structured import BaseStructuredCalibration
from .aperture import GeometricAperture
from .traces import to_ds9_reg as to_ds9_reg_function


class GeometricModel(GeometricAperture):
    def __init__(self, fibid, boxid, start, stop, model):
        super(GeometricModel, self).__init__(fibid, boxid, start, stop)
        self.model = model

    @property
    def valid(self):
        return self.is_valid()

    def is_valid(self):
        if self.model:
            return True
        else:
            return False

    def __getstate__(self):
        state = super(GeometricModel, self).__getstate__()
        # del state['model']

        params = state['model']['params']
        newparams = {}
        for key, val in params.items():
            serial = json_serial_function(val)
            newparams[key] = serial

        state['model']['params'] = newparams
        return state

    def __setstate__(self, state):
        super(GeometricModel, self).__setstate__(state)
        self._set_model(state['model'])

    def _set_model(self, model):
        if model:
            params = {}
            for key, val in model['params'].items():
                params[key] = convert_function(val)
            model['params'] = params

    @property
    def polynomial(self):
        # FIXME: this is a workaround
        return self.model['params']['mean']

    def aper_center(self):
        return self.model['params']['mean']


class ModelMap(BaseStructuredCalibration):

    __tags__ = ['insmode', 'vph']

    def __init__(self, instrument='unknown'):
        super(ModelMap, self).__init__(instrument)
        self.contents = []
        self.boxes_positions = []
        self.global_offset = nppol.Polynomial([0.0])
        self.ref_column = 2000
        self._wcols = None

    def __getstate__(self):
        st = super(ModelMap, self).__getstate__()
        st['contents'] = [t.__getstate__() for t in self.contents]
        st['boxes_positions'] = self.boxes_positions
        st['global_offset'] = self.global_offset.coef
        st['ref_column'] = self.ref_column
        return st

    def __setstate__(self, state):
        super(ModelMap, self).__setstate__(state)
        # self.contents = [GeometricModel(**trace) for trace in state['contents']]
        self.contents = []
        for trace in state['contents']:
            m = GeometricModel.__new__(GeometricModel)
            m.__setstate__(trace)
            self.contents.append(m)
        self.boxes_positions = state.get('boxes_positions', [])
        self.global_offset = nppol.Polynomial(
            state.get('global_offset', [0.0]))
        self.ref_column = state.get('ref_column', 2000)
        self._wcols = None

    def calculate_matrices(self, shape, processes=0):
        if self._wcols is None:
            self._wcols = calc_matrix_cols(self, shape, processes)

    def aper_extract(self, img, processes=0):
        if self._wcols is None:
            self._wcols = calc_matrix_cols(self, img.shape, processes)
        return aper_extract(self, self._wcols, img)

    def to_ds9_reg(self, ds9reg, rawimage=False, numpix=100, fibid_at=0):
        """Transform fiber traces to ds9-region format.

        Parameters
        ----------
        ds9reg : BinaryIO
            Handle to output file name in ds9-region format.
        rawimage : bool
            If True the traces must be generated to be overplotted on
            raw FITS images.
        numpix : int
            Number of abscissae per fiber trace.
        fibid_at : int
            Abscissae where the fibid is shown (default=0 -> not shown).

        """
        return to_ds9_reg_function(self, ds9reg,
                                   rawimage=rawimage,
                                   numpix=numpix,
                                   fibid_at=fibid_at
                                   )


# BUILD MATRICES
