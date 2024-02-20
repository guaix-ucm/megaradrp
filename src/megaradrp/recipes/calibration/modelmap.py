#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

""" Trace model recipe for Megara"""

from __future__ import division, print_function

import multiprocessing as mp

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from numina.core import Result, Parameter
from numina.array import combine
from numina.frame.utils import copy_img
from numina.util.objimport import import_object

from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.products.modelmap import ModelMap
from megaradrp.products.modelmap import GeometricModel
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.modelmap import calc1d_model
from megaradrp.processing.modeldesc import config
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
        reduced = basic_processing_with_combination(
            rinput, flow1, method=combine.median)
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

        # model name for fitting (this could be a parameter
        model_name = 'gaussbox'

        model_map = ModelMap(instrument=obresult.instrument)

        self.logger.debug('update metadata in model')
        model_map.update_metadata(self)
        fp_conf = FocalPlaneConf.from_img(reduced)
        model_map.total_fibers = fp_conf.nfibers
        model_map.missing_fibers = rinput.master_traces.missing_fibers
        model_map.tags = self.extract_tags_from_ref(
            reduced, model_map.tag_names(), base=obresult.labels)
        # model_map.boxes_positions = box_borders
        # model_map.ref_column = cstart
        model_map.update_metadata(self)
        model_map.update_metadata_origin(obresult_meta)
        # Temperature in Celsius with 2 decimals
        model_map.tags['temp'] = round(
            obresult_meta['info'][0]['temp'] - 273.15, 2)

        self.logger.info('perform model fitting')

        tracemap = rinput.master_traces
        data = reduced[0].data

        cols = range(100, 4100, 100)
        # cols = range(100, 200, 100)
        valid_fibers = [(f.fibid, f.boxid)
                        for f in tracemap.contents if f.valid]
        # ncol = tracemap.total_fibers

        nfit = data.shape[1]

        if model_name == 'gaussbox':
            # parameters for this model
            sigma = 1.53
            model_kwargs = {'sigma': sigma}
        else:
            raise ValueError(f'model name {model_name} is undefined')

        objpath = config[model_name]
        model_class = import_object(objpath)
        model_obj = model_class(**model_kwargs)

        # Perform fitting with multiprocessing
        results_get = fit_model(model_obj, data, tracemap, cols,
                                processes=processes)

        self.logger.info('perform model fitting end')

        self.logger.info('interpolate parameters')

        # summarize values
        params_save = model_obj.params_save
        params = model_obj.params_fit
        spline_degrees = model_obj.deg_save

        # vector of columns where we have performed the fit
        g_col = np.asarray(cols)

        for fibid, boxid in valid_fibers:
            # interpolator of the parameters of a given fiber
            interpolators = {name: None for name in params_save}
            # Values of the parameters of a given fiber
            g_vals = {name: [] for name in params}

            # log only 1 in 100 fibers
            dolog = (fibid % 100 == 0)

            if dolog:
                self.logger.debug('compute fibid %d', fibid)

            for calc_col, vals in results_get:
                # Parameters in a given column and fiber fibid
                param_col_fib = vals[fibid]
                for name in params:
                    g_vals[name].append(param_col_fib[name])

            # Fit a UnivariateSpline to each storable parameter
            for name, deg in zip(params_save, spline_degrees):
                interpolators[name] = UnivariateSpline(
                    g_col, g_vals[name], k=deg)

            if self.intermediate_results:
                if dolog:
                    self.logger.debug('creating plots')
                # plot each storable parameter
                for name in params_save:
                    plt.title(f'{name} fib{fibid:03d}')
                    plt.plot(g_col, g_vals[name], 'b*')
                    plt.plot(g_col, interpolators[name](g_col), 'r')
                    plt.savefig(f'fib_{fibid:03d}_{name}.png')
                    plt.close()

                if dolog:
                    self.logger.debug('creating plots end')

            # summary of model for this fiber
            model_fib = {'model_name': model_obj.name, 'params': interpolators}
            # if invalid. missing, model = {}
            gm = GeometricModel(fibid, boxid,
                                start=1, stop=nfit, model=model_fib
                                )

            model_map.contents.append(gm)

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


def calc_parallel(model_desc, data, calc_col, tracemap,
                  nloop=10, average=0):

    if average > 0:
        column = data[:, calc_col -
                      average:calc_col - average + 1].mean(axis=1)
    else:
        column = data[:, calc_col]

    valid_fibers = [f.fibid for f in tracemap.contents if f.valid]
    centers = np.array([f.polynomial(calc_col)
                       for f in tracemap.contents if f.valid])
    # we might need a better approach to logging in multiprocessing
    # https://www.jamesfheath.com/2020/06/logging-in-python-while-multiprocessing.html
    print('computing in column', calc_col)

    # Scale image value
    scale = column.max()
    column_norm = column / scale

    final = calc1d_model(model_desc, column_norm, centers,
                         valid_fibers, calc_col, lateral=2, nloop=nloop)

    # TODO: we may need a function to perform scaling in general
    for idx, params in final.items():
        final[idx]['amplitude'] *= scale

    return calc_col, final


def fit_model(model_desc, data, tracemap, cols, processes=20):

    pool = mp.Pool(processes)

    results = [pool.apply_async(
        calc_parallel,
        args=(model_desc, data, col, tracemap),
        kwds={'nloop': 3, 'average': 2}
    ) for col in cols]

    results_get = [p.get() for p in results]
    return results_get
