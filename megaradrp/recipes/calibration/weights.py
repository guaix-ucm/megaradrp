#
# Copyright 2015-2017 Universidad Complutense de Madrid
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

""" Extraction weights calibration recipes for Megara"""

from __future__ import division

import copy
import types
import math
import multiprocessing as mp
import json as ujson
from tempfile import mkdtemp
import os
import os.path

import six.moves.copyreg as copyreg
import numpy as np
from astropy.modeling import fitting
from scipy.stats import norm
from astropy.modeling.models import custom_model
from numina.core import Product, Requirement
from numina.array import combine

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.types import MasterWeights


copyreg.pickle(
    types.MethodType,
    lambda method: (getattr, (method.im_self, method.im_func.__name__)),
    getattr
)


M_SQRT_2_PI = math.sqrt(2 * math.pi)


class WeightsRecipe(MegaraBaseRecipe):
    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_fiberflat_frame = reqs.MasterFiberFlatFrameRequirement()
    tracemap = reqs.MasterTraceMapRequirement()

    # Products
    master_weights = Product(MasterWeights)

    def __init__(self, size=4096, fibers=623, rows=4112, *args, **kwargs):
        self.SIZE = size
        self.ROWS = rows
        self.FIBERS = fibers
        self.procesos = mp.cpu_count() - 2

        super(WeightsRecipe, self).__init__(version="0.1.0", *args, **kwargs)

    def _add_file_to_tar(self, file_name, tar):
        '''
        :param file_name: <str> Name of the *.fits files
        :param tar: <tarfile> descriptor of the tarfile object
        :return:
        '''
        tar.add(file_name, arcname=os.path.basename(file_name))

    def _check_directory(self, path):
        '''
        :param path: <str> Path where fits files are stored. If exists then will be erased
        :return: None
        '''
        import shutil
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)

    def compress(self, path, tar_name='files'):
        '''
        :param path: path: <str> Path where tar file is stored
        :param tar_name: <str> name of the tar file
        :return: None
        '''
        import tarfile
        import glob

        try:
            os.remove("%s.tar" % tar_name)
        except OSError:
            pass
        tar = tarfile.open("%s.tar" % tar_name, "w")
        files = glob.glob('%s/*.*' % path)
        for file in files:
            self._add_file_to_tar(file, tar)
        tar.close()

    def decompress(self, tar_name='files'):
        '''
        :param tar_name: <str> name of the tar file
        :return: None
        '''
        import tarfile

        tar = tarfile.open("%s.tar" % tar_name, 'r')
        aux = tar.extractall()
        try:
            return tar.getnames()[0].split('/')[0]
        except:
            return ''

    def extract_w(self, img, mlist=[]):
        '''
        :param img: <fits> original fiber flat fits file
        :param mlist: <list> list of csr_matrix
        :return: <ndarray> result of lsqr
        '''
        from scipy.sparse.linalg import lsqr

        result = np.zeros((self.FIBERS, self.SIZE))
        for col in range(self.SIZE):
            wes_csr = mlist[col]
            p = img[:, col]
            x = lsqr(wes_csr, p)
            result[:, col] = x[0]
        return result

    def _load_files_paralell(self, col, path):
        '''
        :param col: <str,int> name of the fits file. It is a counter
        :param path: <str> path where *.npz are
        :return: csr_matrix
        '''
        from scipy.sparse import csr_matrix

        filename = '%s/%s.npz' % (path, col)
        loader = np.load(filename)
        return csr_matrix(
            (loader['data'], loader['indices'], loader['indptr']),
            shape=loader['shape'])

    def load_files_from_directory(self, path, tar_file=None):
        '''
        :param path: <str> path to load *.fits files
        :param tar_file: <str> if it is given, *.fits files are extracted
        :return: list of csr_matrix
        '''

        if tar_file:
            path = self.decompress()

        pool = mp.Pool(processes=self.procesos)
        results = [pool.apply_async(self._load_files_paralell,
                                    args=(ite, path)) for ite in
                   range(self.SIZE)]
        results = [p.get() for p in results]
        return results

    def pixcont(self, i, x0, sig, hpix=0.5):
        '''Integrate a gaussian profile.'''
        z = (i - x0) / sig
        hpixs = hpix / sig
        z2 = z + hpixs
        z1 = z - hpixs
        return norm.cdf(z2) - norm.cdf(z1)

    def g_profile(self, xl, l, s):
        '''A gaussian profile.'''
        z = (xl - l) / s
        return np.exp(-0.5 * z ** 2)

    def fit1d_profile(self, xl, yl, init0, N, nloop=10, S=3):
        """Iterative fitting"""

        init = copy.deepcopy(init0)

        changes_a = np.zeros((N, nloop))
        changes_m = np.zeros((N, nloop))
        changes_s = np.zeros((N, nloop))

        for il in range(nloop):

            values = np.random.permutation(N)

            for val in values:

                m1 = max(0, int(init[val]['mean']) - 6 * S)
                m2 = int(init[val]['mean']) + 6 * S

                y = yl[m1:m2].copy()
                xt = xl[m1:m2]

                for peakid in range(max(0, val - S), min(N, val + S + 1)):
                    if peakid == val:
                        continue

                    y -= gauss_box_model(xt, **init[peakid])

                model = GaussBox(**init[val])
                model.mean.min = model.mean.value - 0.5
                model.mean.max = model.mean.value + 0.5
                # model.mean.fixed = True
                model.stddev.min = 1.0
                model.stddev.max = 2.0
                model.hpix.fixed = True

                fitter = fitting.LevMarLSQFitter()
                model_fitted = fitter(model, xt, y)

                na = model_fitted.amplitude.value
                nm = model_fitted.mean.value
                ns = model_fitted.stddev.value

                changes_a[val, il] = na - init[val]['amplitude']
                changes_m[val, il] = nm - init[val]['mean']
                changes_s[val, il] = ns - init[val]['stddev']

                init[val]['amplitude'] = na
                init[val]['mean'] = nm
                init[val]['stddev'] = ns

        return init, (changes_a, changes_m, changes_s)

    def calc_sparse_matrix(self, final, nrows, cut=1.0e-6, extra=10):
        from scipy.sparse import lil_matrix

        idxs = range(len(final))

        #    g_ampl = np.array([final[i]['amplitude'] for i in idxs])
        g_mean = np.array([final[i]['mean'] for i in idxs])
        g_std = np.array([final[i]['stddev'] for i in idxs])

        # calc w
        begpix = np.ceil(g_mean - 0.5).astype('int')

        steps = np.arange(-extra, extra)
        ref = begpix + steps[:, np.newaxis]

        rr = gauss_box_model(ref, mean=g_mean, stddev=g_std)
        rrb = begpix - extra
        # Filter values below 'cut'
        rr[rr < cut] = 0.0

        # Calc Ws matrix
        block, nfib = rr.shape
        w_init = lil_matrix((nrows, nfib))

        for i in range(nfib):
            w_init[rrb[i]:rrb[i] + block, i] = rr[:, i, np.newaxis]

        # Convert to CSR matrix
        wcol = w_init.tocsr()
        return wcol

    def calc_profile(self, data1, pols, col, sigma, start=0, doplots=False):
        # print 'calc_profile: fitting column', col

        peaks = np.array([pol(col) for pol in pols])

        boxd = data1[:, col]

        centers = peaks[:] - start
        sigs = sigma * np.ones_like(centers)
        scale_sig = 0.25  # For sigma ~= 1.5, the peak is typically 0.25
        ecenters = np.ceil(centers - 0.5).astype('int')

        N = len(centers)
        cmax = boxd.max()
        yl = boxd / cmax  # Normalize to peak
        xl = np.arange(len(yl))

        init_vals = {}
        for i in range(N):
            init_vals[i] = {}
            init_vals[i]['amplitude'] = yl[ecenters[i]] / scale_sig
            # init_vals[i]['mean'] = ecenters[i]
            init_vals[i]['mean'] = centers[i]
            init_vals[i]['stddev'] = sigma

        final, changes = self.fit1d_profile(xl, yl, init_vals, N, nloop=10)

        for i in range(N):
            final[i]['amplitude'] = final[i]['amplitude'] * cmax

        return final

    def run(self, rinput):
        temporary_path = mkdtemp()

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        data2 = img

        pols2 = [t.polynomial for t in rinput.tracemap.contents]

        nrows = data2[0].shape[0]  # 4112
        total_number = data2[0].shape[1]
        cols = range(total_number)  # 4096   # ORIGINAL

        self._check_directory(os.path.join(temporary_path,'chunks'))
        self._check_directory(os.path.join(temporary_path,'json'))

        pool = mp.Pool(processes=self.procesos)
        results = [pool.apply_async(calc_all, args=(ite, data2[0].data, pols2, nrows, temporary_path)) for ite in cols]
        results = [p.get() for p in results]

        self.compress(os.path.join(temporary_path,'chunks'),os.path.join(temporary_path,'master_weights'))
        result = self.create_result(master_weights=os.path.join(temporary_path,'master_weights.tar'))
        # shutil.rmtree(temporary_path)
        return result


# FIXME: GaussBox is duplicated here


def norm_pdf_t(x):
    return np.exp(-0.5 * x * x) / M_SQRT_2_PI


def gauss_box_model_deriv(x, amplitude=1.0, mean=0.0, stddev=1.0, hpix=0.5):
    '''Integrate a gaussian profile.'''

    z = (x - mean) / stddev
    z2 = z + hpix / stddev
    z1 = z - hpix / stddev

    da = norm.cdf(z2) - norm.cdf(z1)

    fp2 = norm_pdf_t(z2)
    fp1 = norm_pdf_t(z1)

    dl = -amplitude / stddev * (fp2 - fp1)
    ds = -amplitude / stddev * (fp2 * z2 - fp1 * z1)
    dd = amplitude / stddev * (fp2 + fp1)

    return (da, dl, ds, dd)


def gauss_box_model(x, amplitude=1.0, mean=0.0, stddev=1.0,
                        hpix=0.5):
    '''Integrate a gaussian profile.'''
    z = (x - mean) / stddev
    m2 = z + hpix / stddev
    m1 = z - hpix / stddev
    return amplitude * (norm.cdf(m2) - norm.cdf(m1))


def pixcont(i, x0, sig, hpix=0.5):
    '''Integrate a gaussian profile.'''
    z = (i - x0) / sig
    hpixs = hpix / sig
    z2 = z + hpixs
    z1 = z - hpixs
    return norm.cdf(z2) - norm.cdf(z1)


def g_profile(xl, l, s):
    '''A gaussian profile.'''
    z = (xl - l) / s
    return np.exp(-0.5 * z ** 2)


def fit1d_profile(xl, yl, init0, N, nloop=10, S=3):
    """Iterative fitting"""

    init = copy.deepcopy(init0)

    changes_a = np.zeros((N, nloop))
    changes_m = np.zeros((N, nloop))
    changes_s = np.zeros((N, nloop))

    for il in range(nloop):

        values = np.random.permutation(N)

        for val in values:

            m1 = max(0, int(init[val]['mean']) - 6 * S)
            m2 = int(init[val]['mean']) + 6 * S

            y = yl[m1:m2].copy()
            xt = xl[m1:m2]

            for peakid in range(max(0, val - S), min(N, val + S + 1)):
                if peakid == val:
                    continue

                y -= gauss_box_model(xt, **init[peakid])

            model = GaussBox(**init[val])
            model.mean.min = model.mean.value - 0.5
            model.mean.max = model.mean.value + 0.5
            # model.mean.fixed = True
            model.stddev.min = 1.0
            model.stddev.max = 2.0
            model.hpix.fixed = True

            fitter = fitting.LevMarLSQFitter()
            model_fitted = fitter(model, xt, y)

            na = model_fitted.amplitude.value
            nm = model_fitted.mean.value
            ns = model_fitted.stddev.value

            changes_a[val, il] = na - init[val]['amplitude']
            changes_m[val, il] = nm - init[val]['mean']
            changes_s[val, il] = ns - init[val]['stddev']

            init[val]['amplitude'] = na
            init[val]['mean'] = nm
            init[val]['stddev'] = ns

    return init, (changes_a, changes_m, changes_s)


def calc_sparse_matrix(final, nrows, cut=1.0e-6, extra=10):
    from scipy.sparse import lil_matrix

    idxs = range(len(final))

    #    g_ampl = np.array([final[i]['amplitude'] for i in idxs])
    g_mean = np.array([final[i]['mean'] for i in idxs])
    g_std = np.array([final[i]['stddev'] for i in idxs])

    # calc w
    begpix = np.ceil(g_mean - 0.5).astype('int')

    steps = np.arange(-extra, extra)
    ref = begpix + steps[:, np.newaxis]

    rr = gauss_box_model(ref, mean=g_mean, stddev=g_std)
    rrb = begpix - extra
    # Filter values below 'cut'
    rr[rr < cut] = 0.0

    # Calc Ws matrix
    block, nfib = rr.shape
    w_init = lil_matrix((nrows, nfib))

    for i in range(nfib):
        w_init[rrb[i]:rrb[i] + block, i] = rr[:, i, np.newaxis]

    # Convert to CSR matrix
    wcol = w_init.tocsr()
    return wcol


def calc_profile(data1, pols, col, sigma, start=0, doplots=False):
    # print 'calc_profile: fitting column', col

    peaks = np.array([pol(col) for pol in pols])

    boxd = data1[:, col]

    centers = peaks[:] - start
    sigs = sigma * np.ones_like(centers)
    scale_sig = 0.25  # For sigma ~= 1.5, the peak is typically 0.25
    ecenters = np.ceil(centers - 0.5).astype('int')

    N = len(centers)
    cmax = boxd.max()
    yl = boxd / cmax  # Normalize to peak
    xl = np.arange(len(yl))

    init_vals = {}
    for i in range(N):
        init_vals[i] = {}
        init_vals[i]['amplitude'] = yl[ecenters[i]] / scale_sig
        # init_vals[i]['mean'] = ecenters[i]
        init_vals[i]['mean'] = centers[i]
        init_vals[i]['stddev'] = sigma

    final, changes = fit1d_profile(xl, yl, init_vals, N, nloop=10)

    for i in range(N):
        final[i]['amplitude'] = final[i]['amplitude'] * cmax

    return final


def calc_all(col, data2, pols2, nrows, temporary_path):
    '''
    Poner bien los directorios
    :param col:
    :return:
    '''
    prefix = os.path.join(temporary_path,'json')
    sigma = 1.5  # Typical value for MEGARA
    fname = os.path.join(prefix, '%d.json' % (col,))

    final = calc_profile(data2, pols2, col, sigma, start=0)

    with open(fname, 'w') as outfile:
        ujson.dump(final, outfile)

    wm = calc_sparse_matrix(final, nrows, cut=1.0e-6, extra=10)

    prefixw = os.path.join(temporary_path,'chunks')
    jsonname = os.path.join(prefixw, '%d' % (col,))

    np.savez(jsonname, data=wm.data, indices=wm.indices, indptr=wm.indptr,
             shape=wm.shape)

    return final, col


GaussBox = custom_model(gauss_box_model,func_fit_deriv=gauss_box_model_deriv)
