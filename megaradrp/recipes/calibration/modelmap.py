#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

""" Trace model recipe for Megara"""

from __future__ import division, print_function

import math
import bisect
import multiprocessing as mp

import numpy as np
from astropy.modeling import fitting
from numina.core import Product, Requirement, Parameter
from numina.array import combine
from numina.modeling.gaussbox import GaussBox, gauss_box_model
from scipy.interpolate import UnivariateSpline

from megaradrp.products.modelmap import ModelMap
from megaradrp.products.modelmap import GeometricModel
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.types import ProcessedImage, ProcessedRSS
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs


class ModelMapRecipe(MegaraBaseRecipe):
    # Requirements
    master_bpm = reqs.MasterBPMRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    # FIXME: this is not really necessary, it can be computed
    # from the data
    master_traces = reqs.MasterTraceMapRequirement()
    processes = Parameter(0, 'Number of processes used for fitting')
    debug_plot = Parameter(0, 'Save intermediate tracing plots')
    # Products
    reduced_image = Product(ProcessedImage)
    reduced_rss = Product(ProcessedRSS)
    master_model = Product(ModelMap)

    def run(self, rinput):

        self.logger.info('starting model map recipe')

        obresult = rinput.obresult
        obresult_meta = obresult.metadata_with(self.datamodel)
        if rinput.processes == 0:
            have = mp.cpu_count()
            if have >= 4:
                processes = mp.cpu_count() - 2
            else:
                processes = 1
        else:
            processes = rinput.processes
        self.logger.debug('using %d processes', processes)

        self.logger.info('start basic reduction')
        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        reduced = basic_processing_with_combination(rinput, flow1, method=combine.median)
        self.set_base_headers(reduced[0].header)
        self.logger.info('end basic reduction')

        self.save_intermediate_img(reduced, 'reduced_image.fits')

        self.logger.debug('create result model')

        # insconf = obresult.configuration

        # boxes = insconf.get('pseudoslit.boxes', **obresult.tags)

        # box_borders0, cstart0 = self.obtain_boxes(insconf, obresult.tags)

        # box_borders, cstart = self.refine_boxes_from_image(reduced, box_borders0, cstart0)

        # self.logger.debug("original boxes: %s", box_borders0)
        # self.logger.debug("refined boxes: %s", box_borders)

        model_map = ModelMap(instrument=obresult.instrument)

        self.logger.debug('update metadata in model')
        model_map.update_metadata(self)
        fiberconf = self.datamodel.get_fiberconf(reduced)
        model_map.total_fibers = fiberconf.nfibers
        model_map.missing_fibers = rinput.master_traces.missing_fibers
        model_map.tags = obresult.tags
        # model_map.boxes_positions = box_borders
        # model_map.ref_column = cstart
        model_map.update_metadata(self)
        model_map.update_metadata_origin(obresult_meta)
        # Temperature in Celsius with 2 decimals
        model_map.tags['temp'] = round(obresult_meta['info'][0]['temp'] - 273.15, 2)

        self.logger.info('perform model fitting')

        tracemap = rinput.master_traces
        data = reduced[0].data
        sigma = 1.53
        cols = range(100, 4100, 100)
        # cols = range(100, 200, 100)
        valid = [f.fibid for f in tracemap.contents if f.valid]
        ncol = tracemap.total_fibers
        nrow = data.shape[0]
        nfit = data.shape[1]
        results_get = fit_model(data, tracemap, valid, nrow, ncol, sigma, cols, processes)

        self.logger.info('perform model fitting end')

        self.logger.info('interpolate parameters')

        array_mean = np.zeros((ncol, nfit))
        array_std = np.zeros((ncol, nfit))

        xcol = np.arange(nfit)

        mean_splines = {}
        std_splines = {}

        for fibid in valid:
            g_amp = []
            g_std = []
            g_mean = []
            g_col = []

            for calc_col, vals in results_get:
                param = vals[fibid]
                g_col.append(calc_col)
                g_std.append(param['stddev'])
                g_mean.append(param['mean'])
                g_amp.append(param['amplitude'])

            interpol_std = UnivariateSpline(g_col, g_std, k=5)
            interpol_mean = UnivariateSpline(g_col, g_mean, k=3)

            mean_splines[fibid] = interpol_mean
            std_splines[fibid] = interpol_std

            row = fibid - 1
            array_mean[row] = interpol_mean(xcol)
            array_std[row] = interpol_std(xcol)

            params = {'stddev': interpol_std, 'mean': interpol_mean}
            model = {'model_name': 'gaussbox', 'params': params}
            # if invalid. missing, model = {}
            m = GeometricModel(
                fibid,
                boxid=1, # FIXME: not counting this
                start=1,
                stop=nfit,
                model=model
            )

            model_map.contents.append(m)
            if False:
                import matplotlib.pyplot as plt
                plt.title('std %s' % fibid)
                plt.plot(g_col, g_std, 'b*')
                plt.plot(g_col, interpol_std(g_col), 'r')
                plt.show()

                plt.title('position %s' % fibid)
                plt.plot(g_col, g_mean, 'b*')
                plt.plot(g_col, interpol_mean(g_col), 'r')
                plt.show()

        self.logger.info('interpolate parameters end')

        # perform extraction with our own calibration
        self.logger.info('perform extraction with computed calibration')
        calibrator_aper = ApertureExtractor(model_map, self.datamodel)
        reduced_rss = calibrator_aper(reduced)

        self.logger.info('ending model map recipe')
        result = self.create_result(reduced_image=reduced,
                                    reduced_rss=reduced_rss,
                                    master_model=model_map)

        return result


def initial_base(boxd1d, ecenters, npix=5):
    nfib = len(ecenters)

    offset = npix // 2
    yfit = np.empty((npix, nfib))
    xfit = np.arange(npix) - offset

    for i in range(npix):
        yfit[i, :] = boxd1d[ecenters + (i - offset)]

    coeff = np.polyfit(xfit, np.log(yfit), deg=2)
    c, b, a = coeff

    sig2 = -1 / (2 * c)
    mu = b / sig2
    ampl = np.exp(a - mu ** 2)

    return ampl, mu, sig2


def calc1d_N(boxd1d, valid, centers1d, sigma, lateral=0, reject=3, nloop=1, fixed_centers=False):
    # initial(boxd1d, valid, centers1d)

    from astropy.modeling.functional_models import Const1D

    ecenters = np.ceil(centers1d - 0.5).astype('int')
    nfib = len(centers1d)
    # FIXME, hardcoded
    xl = np.arange(4112.0)
    ampl, mu, sig2 = initial_base(boxd1d, ecenters, npix=5)
    sig_calc = np.sqrt(sig2)
    yl = boxd1d
    init_simple = True
    init = {}
    for i in range(nfib):
        fibid = valid[i]
        init[fibid] = {}
        if init_simple:
            init[fibid]['amplitude'] = boxd1d[ecenters[i]] / 0.25
            init[fibid]['stddev'] = sigma
            init[fibid]['mean'] = centers1d[i]
        else:
            init[fibid]['amplitude'] = ampl[i]
            init[fibid]['mean'] = centers1d[i]
            init[fibid]['stddev'] = sig_calc[i]

    # plt.plot(sig_calc / sigma)
    # plt.show()

    # total fits is 2 * lateral + 1
    total = reject + lateral

    for il in range(nloop):

        # print 'loop', il, datetime.datetime.now()

        # permutation of the valid fibers
        values = np.random.permutation(valid)

        for val in values:
            #

            m1 = max(0, int(math.ceil(init[val]['mean'] - 0.5)) - 6 * total)
            m2 = int(math.ceil(init[val]['mean'] - 0.5)) + 6 * total

            y = yl[m1:m2].copy()
            xt = xl[m1:m2]

            pos = bisect.bisect(valid, val)
            candidates = valid[max(pos - total - 1, 0): pos + total]
            dis_f = [abs(c - val) for c in candidates]
            # for p in range(max(0, val-S), min(nfib, val+S+1)):
            fit_ids = []
            for c, df in zip(candidates, dis_f):
                if df > lateral:
                    # remove contribution
                    y -= gauss_box_model(xt, **init[c])
                else:
                    fit_ids.append(c)

            num_of_fibers = len(fit_ids)
            # print('lateral', lateral, val, fit_ids)
            # print('num', num_of_fibers)

            Model = Const1D
            for _ in range(num_of_fibers):
                Model = Model + GaussBox

            # Extract values to create initials
            sub_init = {}
            sub_init["amplitude_0"] = 0.0
            for idx, fib_id in enumerate(fit_ids, 1):
                key = "amplitude_%d" % idx
                sub_init[key] = init[fib_id]['amplitude']
                key = "mean_%d" % idx
                sub_init[key] = init[fib_id]['mean']
                key = "stddev_%d" % idx
                sub_init[key] = init[fib_id]['stddev']

            # Model = sum_of_gaussian_factory(num_of_fibers)
            model = Model(**sub_init)

            # Set fit conditions
            model.amplitude_0.fixed = True
            for idx, fib_id in enumerate(fit_ids, 1):
                key = "mean_%d" % idx
                m_mean = getattr(model, key)
                if fixed_centers:
                    m_mean.fixed = True
                else:
                    m_mean.min = init[fib_id]['mean'] - 0.5
                    m_mean.max = init[fib_id]['mean'] + 0.5

                key = "stddev_%d" % idx
                m_stddev = getattr(model, key)
                m_stddev.min = init[fib_id]['stddev'] - 0.5
                m_stddev.max = init[fib_id]['stddev'] + 0.5

                key = "hpix_%d" % idx
                m_hpix = getattr(model, key)
                m_hpix.fixed = True

            fitter = fitting.LevMarLSQFitter()
            model_fitted = fitter(model, xt, y)

            # Extract values
            # This extracts only val
            # val_idx = fit_ids.index(val) + 1
            # na = getattr(model_fitted, 'amplitude_%d' % val_idx).value
            # nm = getattr(model_fitted, 'mean_%d' % val_idx).value
            # ns = getattr(model_fitted, 'stddev_%d' % val_idx).value

            # init[val]['amplitude'] = na
            # init[val]['mean'] = nm
            # init[val]['stddev'] = ns

            for val_idx, fibid in enumerate(fit_ids, 1):
                # This will overwrite fits with fits performed later
                # if fibid == val:
                na = getattr(model_fitted, 'amplitude_%d' % val_idx).value
                nm = getattr(model_fitted, 'mean_%d' % val_idx).value
                ns = getattr(model_fitted, 'stddev_%d' % val_idx).value

                init[fibid]['amplitude'] = na
                init[fibid]['mean'] = nm
                init[fibid]['stddev'] = ns
                # break

    return init


def calc2_base(column, centers, valid, sigma, nloop=10, do_plot=False):
    scale = column.max()
    column_norm = column / scale
    final = calc1d_N(column_norm, valid, centers, sigma, lateral=2, nloop=nloop)

    for idx, params in final.items():
        final[idx]['amplitude'] *= scale

    if do_plot:
        pass # calc2_plot(final, centers, valid, sigma)

    return final


def calc_para_2(data, calc_col, tracemap, valid, nrow, ncol, sigma,
                nloop=10, clip=1.0e-6, extra=10, average=0, calc_init=True):
    if average > 0:
        yl = data[:, calc_col - average:calc_col - average + 1].mean(axis=1)
    else:
        yl = data[:, calc_col]

    centers = np.array([f.polynomial(calc_col) for f in tracemap.contents if f.valid])
    print ('calc_col is', calc_col)

    final = calc2_base(yl, centers, valid, sigma, nloop=nloop, do_plot=False)

    return calc_col, final


def fit_model(data, tracemap, valid, nrow, ncol, sigma, cols, processes=20):

    pool = mp.Pool(processes)

    results = [pool.apply_async(
        calc_para_2,
        args=(data, col, tracemap, valid, nrow, ncol, sigma),
        kwds={'nloop': 3, 'average': 2, 'calc_init': True}
    ) for col in cols]

    results_get = [p.get() for p in results]
    return results_get
