#
# Copyright 2017-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""

import numpy
import numpy.polynomial.polynomial as nppol

from numina.util.convertfunc import json_serial_function, convert_function

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
                params[key] =  convert_function(val)
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
        self.global_offset = nppol.Polynomial(state.get('global_offset', [0.0]))
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
def calc_matrix(g_mean, g_std, valid, wshape, clip=1.0e-6, extra=10):
    from scipy.sparse import lil_matrix
    from numina.modeling.gaussbox import gauss_box_model
    # calc w
    # this is valid for 1 column (there are more broadcasts for N)

    begpix = numpy.ceil(g_mean - 0.5).astype('int')

    steps = numpy.arange(-extra, extra)
    # ref is a bidimensional matrix +-10 pixels around the trace
    ref = begpix + steps[:, numpy.newaxis]

    rr = gauss_box_model(ref, mean=g_mean, stddev=g_std)
    rrb = begpix - extra

    # This was sending warnings. Is there a NaN somewhere?
    with numpy.errstate(invalid='ignore'):
        rr[rr < clip] = 0.0

    # Calc Ws matrix
    block, valid_nfib = rr.shape

    w_init = lil_matrix(wshape)
    for fibid in valid:
        idx = fibid - 1
        w_init[rrb[idx]:rrb[idx] + block, idx] = rr[:, idx, numpy.newaxis]

    wcol = w_init.tocsr()
    return wcol


def calc_matrix_alt(g_mean, g_std, col, valid, wshape, clip=1.0e-6, extra=10):

    # For parallel processing
    wcol = calc_matrix(g_mean, g_std, valid, wshape, clip=clip, extra=extra)

    return col, wcol


def calc_matrix_cols(model_map, datashape, processes=0):

    dnrow, dncol = datashape
    nfibs = model_map.total_fibers

    shape = (nfibs, dncol)
    wshape = (dnrow, nfibs)

    array_mean = numpy.zeros(shape)
    array_std = numpy.zeros(shape)

    xcol = numpy.arange(dncol)
    for pesos in model_map.contents:
        if pesos.valid:
            row = pesos.fibid - 1
            params = pesos.model['params']
            array_mean[row] = params['mean'](xcol)
            array_std[row] = params['stddev'](xcol)

    # shift array_mean according to global_offset
    mask = numpy.zeros((model_map.total_fibers,), dtype='bool')
    valid_r = [(f.fibid - 1) for f in model_map.contents if f.valid]
    mask[valid_r] = True

    mean_at_ref = array_mean[mask, model_map.ref_column]
    offset = model_map.global_offset(mean_at_ref)
    array_mean[mask, :] = array_mean[mask, :] + offset[:, numpy.newaxis]

    wcols = {}
    valid = [f.fibid for f in model_map.contents if f.valid]
    if processes < 2:
        for col in xcol:
            result = calc_matrix(array_mean[:, col], array_std[:, col],
                                 valid, wshape, clip=1e-6, extra=10)
            wcols[col] = result
    else:
        import multiprocessing as mp
        pool = mp.Pool(processes=processes)

        results = [pool.apply_async(
            calc_matrix_alt,
            args=(array_mean[:, col], array_std[:, col], col, valid, wshape),
            kwds={'clip': 1e-6, 'extra': 10}
        ) for col in xcol]

        for p in results:
            col, wcol = p.get()
            wcols[col] = wcol

    return wcols


def aper_extract(model_map, wcols, img):
    from scipy.sparse.linalg import lsqr

    n0 = model_map.total_fibers
    n1 = img.shape[1]
    rss = numpy.zeros((n0, n1))
    for key, val in wcols.items():
        yl = img[:, key]
        res = lsqr(val, yl)
        rss[:, key] = res[0]

    return rss