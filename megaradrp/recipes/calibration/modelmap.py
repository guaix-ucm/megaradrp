#
# Copyright 2015-2020 Universidad Complutense de Madrid
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
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from astropy.modeling import fitting
from astropy.modeling.functional_models import Const1D
from numina.core import Result, Requirement, Parameter
from numina.array import combine
from numina.modeling.gaussbox import GaussBox, gauss_box_model
from numina.frame.utils import copy_img

from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.products.modelmap import ModelMap
from megaradrp.products.modelmap import GeometricModel
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.ntypes import ProcessedImage, ProcessedRSS
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs


class ModelMapRecipe(MegaraBaseRecipe):
    """Provides fiber profile information from continuum flat images.

    This recipe process a set of continuum flat images obtained in
    *Trace Map* mode and returns the fiber profile information required
    to perform advanced fiber extraction in other recipes. The recipe also
    returns the result of processing the input images upto dark correction.

    See Also
    --------
    megaradrp.products.modelmap.ModelMap: description of ModelMap product
    megaradrp.recipes.calibration.tracemap.TraceMapRecipe: description of TraceMap recipe
    megaradrp.instrument.configs: instrument configuration

    Notes
    -----
    Images provided in `obresult` are trimmed and corrected from overscan,
    bad pixel mask (if `master_bpm` is not None), bias and dark current
    (if `master_dark` is not None).
    Images thus corrected are the stacked using the median.

    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result and will be used for
    fiting the profiles of the fibers.

    The approximate central position of the fibers is obtained from
    `master_traces`. Then, for each 100 columns of the reduced image, a vertical
    cut in the image is fitted to a sum of fiber profiles, being the profile
    a gaussian convolved with a square.

    The fits are made in parallel, being the number of processes controlled
    by the parameter `processes`, with the default value of 0 meaning to use
    the number of cores minus 2 if the number of cores is greater or equal to 4,
    one process otherwise.

    After the columns are fitted, the profiles (central position and sigma)
    are interpolated to all columns using splines. The coefficientes of the
    splines are stored in the final `master_traces` object.

    The recipe returns also the RSS obtained by applying advanced extraction
    to `reduced_image`. As an intermediate result, the recipe proceduces DS9
    reg files with the position of the center of the profiles, that can be
    used with raw and reduced images.

    """

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
    # Results
    reduced_image = Result(ProcessedImage)
    reduced_rss = Result(ProcessedRSS)
    master_model = Result(ModelMap)

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

        # boxes = insconf.get_property('pseudoslit.boxes')
        # values = insconf.get_property('pseudoslit.boxes_positions')
        # cstart0 = values['ref_column']
        # box_borders0 = values['positions']

        # box_borders, cstart = self.refine_boxes_from_image(reduced, box_borders0, cstart0)

        # self.logger.debug("original boxes: %s", box_borders0)
        # self.logger.debug("refined boxes: %s", box_borders)

        model_map = ModelMap(instrument=obresult.instrument)

        self.logger.debug('update metadata in model')
        model_map.update_metadata(self)
        fp_conf = FocalPlaneConf.from_img(reduced)
        model_map.total_fibers = fp_conf.nfibers
        model_map.missing_fibers = rinput.master_traces.missing_fibers
        model_map.tags = self.extract_tags_from_ref(reduced, model_map.tag_names(), base=obresult.labels)
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

            dolog = (fibid % 100 == 0)

            if dolog:
                self.logger.debug('compute fibid %d', fibid)

            for calc_col, vals in results_get:
                param = vals[fibid]
                g_col.append(calc_col)
                g_std.append(param['stddev'])
                g_mean.append(param['mean'])
                g_amp.append(param['amplitude'])

            interpol_std = UnivariateSpline(g_col, g_std, k=5)
            interpol_mean = UnivariateSpline(g_col, g_mean, k=3)

            if self.intermediate_results:
                if dolog:
                    self.logger.debug('saving plots')
                plt.title('std fib{:03d}'.format(fibid))
                plt.plot(g_col, g_std, 'b*')
                plt.plot(g_col, interpol_std(g_col), 'r')
                plt.savefig('fib_{:03d}_std.png'.format(fibid))
                plt.close()
                plt.title('mean fin{:03d}'.format(fibid))
                plt.plot(g_col, g_mean, 'b*')
                plt.plot(g_col, interpol_mean(g_col), 'r')
                plt.savefig('fib_{:03d}_mean.png'.format(fibid))
                plt.close()
                if dolog:
                    self.logger.debug('saving plots end')

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

        self.logger.info('interpolate parameters end')

        # perform extraction with our own calibration
        self.logger.info('perform extraction with computed calibration')
        calibrator_aper = ApertureExtractor(model_map, self.datamodel)
        reduced_copy = copy_img(reduced)
        reduced_rss = calibrator_aper(reduced_copy)

        if self.intermediate_results:
            with open('ds9.reg', 'w') as ds9reg:
                model_map.to_ds9_reg(ds9reg, rawimage=False,
                                     numpix=100, fibid_at=2048)

            with open('ds9_raw.reg', 'w') as ds9reg:
                model_map.to_ds9_reg(ds9reg, rawimage=True,
                                     numpix=100, fibid_at=2048)

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


def calc1d_N(boxd1d, valid, centers1d, sigma, lateral=0, reject=3, nloop=1,
             fixed_centers=False, init_simple=False):

    ecenters = np.ceil(centers1d - 0.5).astype('int')
    nfib = len(centers1d)

    # FIXME, hardcoded
    xl = np.arange(4112.0)

    if init_simple:
        ampl, mu, sig2 = initial_base(boxd1d, ecenters, npix=5)
        sig_calc = np.sqrt(sig2)
    else:
        sig_calc = sigma * np.ones((nfib,))
        ampl = [boxd1d[ecenters[i]] / 0.25 for i in range(nfib)]

    init = {}
    for i in range(nfib):
        fibid = valid[i]
        init[fibid] = {}
        init[fibid]['mean'] = centers1d[i]
        init[fibid]['stddev'] = sig_calc[i]
        init[fibid]['amplitude'] = ampl[i]

    # plt.plot(sig_calc / sigma)
    # plt.show()

    # total fits is 2 * lateral + 1
    total = reject + lateral
    yl = boxd1d

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

            # num_of_fibers = len(fit_ids)
            # print('lateral', lateral, val, fit_ids)
            # print('num', num_of_fibers)

            # Extract values to create initials
            model = Const1D(amplitude=0.0)
            for fib_id in fit_ids:
                newg = GaussBox(
                    amplitude=init[fib_id]['amplitude'],
                    mean=init[fib_id]['mean'],
                    stddev=init[fib_id]['stddev']
                )
                model = model + newg

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
                nloop=10, average=0, calc_init=True):
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
        kwds={'nloop': 3, 'average': 2}
    ) for col in cols]

    results_get = [p.get() for p in results]
    return results_get
