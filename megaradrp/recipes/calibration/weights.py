#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

""" Extraction weights calibration recipes for Megara"""

from __future__ import division, print_function

import copy
import types
import math
import multiprocessing as mp
import json
from tempfile import mkdtemp
import os
import os.path

import six.moves.copyreg as copyreg
import numpy as np
from astropy.modeling import fitting
from scipy.stats import norm
from numina.core import Product, Requirement
from numina.array import combine
from numina.modeling.gaussbox import GaussBox, gauss_box_model
from numina.user.helpers import make_sure_path_exists

from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.types import WeightsMap


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
    tracemap = reqs.MasterTraceMapRequirement()

    # Products
    master_weights = Product(WeightsMap)

    def __init__(self, size=4096, fibers=623, rows=4112, *args, **kwargs):
        self.SIZE = size
        self.ROWS = rows
        self.FIBERS = fibers
        self.procesos = 2 # mp.cpu_count() - 2
        self.work_dir = None
        self.clean_work_dir = False

        super(WeightsRecipe, self).__init__(version="0.1.0", *args, **kwargs)

    def _add_file_to_tar(self, file_name, tar):
        """
        :param file_name: <str> Name of the *.fits files
        :param tar: <tarfile> descriptor of the tarfile object
        :return:
        """
        tar.add(file_name, arcname=os.path.basename(file_name))

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

    def extract_w(self, img, mlist):
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

    def fit1d_profile(self, xl, yl, init0, N, nloop=10, scale=3):
        return fit1d_profile(xl, yl, init0, N, nloop, scale)

    def calc_sparse_matrix(self, final, nrows, cut=1.0e-6, extra=10):
        return calc_sparse_matrix(final, nrows, cut, extra)

    def calc_profile(self, data1, pols, col, sigma, start=0, doplots=False):
        return calc_profile(data1, pols, col, sigma, start, doplots)

    def tear_down(self):
        if self.clean_work_dir and self.work_dir:
            pass

    def run(self, rinput):

        self.logger.info('starting weights map recipe')

        # We need directories for intermediate results storage
        self.work_dir = self.runinfo['work_dir']
        if self.work_dir is None:
            # FIXME: in this case I should delete it
            # tear_down method
            self.work_dir = mkdtemp()
            self.clean_work_dir = True

        # ensure the directories exist
        chunks_dir = os.path.join(self.work_dir, 'chunks')
        make_sure_path_exists(chunks_dir)
        json_dir = os.path.join(self.work_dir, 'json')
        make_sure_path_exists(json_dir)
        self.logger.debug("chunks dir: %s", chunks_dir)
        self.logger.debug("json dir: %s", json_dir)
        result_path1 = os.path.join(self.runinfo['work_dir'], 'master_weights')
        result_path2 = os.path.join(self.runinfo['work_dir'], 'master_weights.tar')

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        self.logger.info('start weigth extraction')
        self.logger.debug('number of processes: %d', self.procesos)
        pols2 = [t.polynomial for t in rinput.tracemap.contents]

        nrows = img[0].shape[0]  # 4112
        total_number = img[0].shape[1]
        cols = range(total_number)  # 4096   # ORIGINAL
        cols = range(100, 4096, 400)

        pool = mp.Pool(processes=self.procesos)
        # results = [pool.apply_async(calc_all, args=(ite, img[0].data, pols2, nrows, temporary_path)) for ite in cols]
        # results_get = [p.get() for p in results]
        # Not async
        results_get = [calc_all(ite, img[0].data, pols2, nrows, json_dir, chunks_dir) for ite in cols]

        self.compress(chunks_dir, result_path1)
        self.logger.info('ending weights map recipe')
        result = self.create_result(master_weights=result_path2)

        return result


def fit1d_profile(xl, yl, init0, N, nloop=10, scale=3):
    """Iterative fitting"""

    init = copy.deepcopy(init0)

    changes_a = np.zeros((N, nloop))
    changes_m = np.zeros((N, nloop))
    changes_s = np.zeros((N, nloop))

    for il in range(nloop):

        values = np.random.permutation(N)

        for val in values:

            m1 = max(0, int(init[val]['mean']) - 6 * scale)
            m2 = int(init[val]['mean']) + 6 * scale

            y = yl[m1:m2].copy()
            xt = xl[m1:m2]

            for peakid in range(max(0, val - scale), min(N, val + scale + 1)):
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

    print('calc_sparse_matrix, nrows:', nrows)
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

    print("p5 init", w_init.shape, rr.shape)
    for i in range(nfib):
        s = w_init[rrb[i]:rrb[i] + block, i].shape
        if s[0] < 20:
            print("w", w_init[rrb[i]:rrb[i] + block, i].shape)
            print("fib, fibid", i, i + 1)
        #w_init[rrb[i]:rrb[i] + block, i] = rr[:, i, np.newaxis]
    print("p6")
    # Convert to CSR matrix
    wcol = w_init.tocsr()
    print("p7")
    return wcol


def calc_profile(data1, pols, col, sigma, start=0, doplots=False):
    print('calc_profile: fitting column', col)

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

    print('calc_profile: init')
    init_vals = {}
    for i in range(N):
        init_vals[i] = {}
        init_vals[i]['amplitude'] = yl[ecenters[i]] / scale_sig
        # init_vals[i]['mean'] = ecenters[i]
        init_vals[i]['mean'] = centers[i]
        init_vals[i]['stddev'] = sigma

    print('calc_profile: do fitting')
    final, changes = fit1d_profile(xl, yl, init_vals, N, nloop=10)

    # Rescale
    for i in range(N):
        final[i]['amplitude'] = final[i]['amplitude'] * cmax
    print('final calc_profile', col)
    return final


def calc_all(col, data2, pols2, nrows, json_path, chunk_path):
    print("json path", col, json_path)
    print("chunk path", col, chunk_path)
    # prefix = os.path.join(temporary_path,'json')
    sigma = 1.5  # Typical value for MEGARA
    fname = os.path.join(json_path, '%d.json' % (col,))
    print("temp file", col, fname)

    print("calc profile", col)
    final = calc_profile(data2, pols2, col, sigma, start=0)

    with open(fname, 'w') as outfile:
        print("dump profile", col)
        json.dump(final, outfile)

    print("calc sparse matrix", col)
    wm = calc_sparse_matrix(final, nrows, cut=1.0e-6, extra=10)

    # prefixw = os.path.join(temporary_path,'chunks')
    npname = os.path.join(chunk_path, '%d' % (col,))

    print("save sparse matrix", col)
    np.savez(npname, data=wm.data, indices=wm.indices, indptr=wm.indptr,
             shape=wm.shape)

    return final, col
