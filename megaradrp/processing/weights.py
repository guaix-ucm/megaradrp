#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging
import multiprocessing as mp
from astropy.io import fits

import numpy as np
from numina.processing import Corrector

_logger = logging.getLogger('numina.processing')


class WeightsCorrector(Corrector):
    '''A Node that corrects from twilight.'''

    def __init__(self, master_weights, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):
        #tagger = TagFits('NUM-MFF', 'MEGARA master_weights correction')

        super(WeightsCorrector, self).__init__(datamodel=datamodel, dtype=dtype)

        self.master_weights = master_weights
        self.processes = mp.cpu_count() - 2
        self.SIZE = 4096

    def decompress(self):
        '''
        :param tar_name: <str> name of the tar file
        :return: None
        '''

        name = self.master_weights.fileobj.name.split('.tar')[0]

        aux = self.master_weights.extractall(name + '/')
        return name

    def run(self, img):
        '''
        :param img: <numpy.array> reduced image
        :param tar_file: <pointer> tar file to be extracted
        :return: list of extracted weights
        '''
        img = img[0].data
        _logger.debug('correct from weights in image ')

        path = self.decompress()

        _logger.info('decompress done')
        _logger.info('Starting: _load_files_paralell')
        pool = mp.Pool(processes=self.processes)
        results = [pool.apply_async(_load_files_paralell,
                                    args=(ite, path)) for ite in
                   range(self.SIZE)]
        results = [p.get() for p in results]
        # return results

        _logger.info('Starting: extract_w_paralell')

        pool2 = mp.Pool(processes=self.processes)
        extracted_w = [pool2.apply_async(extract_w_paralell,
                                         args=(img[:, ite], results[ite])) for
                       ite in range(self.SIZE)]
        extracted_w = [p.get() for p in extracted_w]

        _logger.info('extracted')

        hdu = fits.PrimaryHDU(np.array(extracted_w).T)

        return fits.HDUList([hdu])


def extract_w_paralell(img, mlist):
    '''
    :param img: <fits>
    :param mlist: <list> one element of the csr_matrix
    :return: <ndarray> result of lsqr
    '''
    from scipy.sparse.linalg import lsqr
    x = lsqr(mlist, img)
    return x[0]


def _load_files_paralell(col, path):
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
