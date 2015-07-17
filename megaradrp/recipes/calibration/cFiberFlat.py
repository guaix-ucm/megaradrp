from __future__ import division, print_function

import logging
import os

import numpy
from astropy.io import fits

from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.array.combine import median as c_median
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from megaradrp.processing import OverscanCorrector, TrimImage
from numina.core.oresult import ObservationResult
from numina.core.dataframe import DataFrame

from megaradrp.products import MasterFiberFlat
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement

from megaradrp.trace.traces import init_traces
from megaradrp.core import apextract_tracemap

from megaradrp.recipes.cBase import BaseRecipe

_logger = logging.getLogger('numina.recipes.megara')


class FiberFlatRecipe(BaseRecipe):
    '''Process FIBER_FLAT images and create MASTER_FIBER_FLAT.'''

    # # Requirements
    # master_bias = MasterBiasRequirement()
    # obresult = ObservationResultRequirement()
    # # Products
    # fiberflat_frame = Product(MasterFiberFlat)
    # fiberflat_rss = Product(MasterFiberFlat)
    # traces = Product(TraceMap)

    def __init__(self, obsrun , workenv, requirement_file):
        super(FiberFlatRecipe,self).__init__(workenv, obsrun, requirement_file)

        # TODO controlar que existe el fichero. Si no, tirar exception
        file = os.path.abspath(os.path.join(self.workenv['workdir'], requirement_file['requirements']['master_bias']))

        mbias = fits.open(file)[0].data.copy()

        self.flow = [OverscanCorrector(), TrimImage(), BiasCorrector(mbias)]

        self.result = []

        #TODO Poner esto mejor y hacerlo con los requirement
        lista = []
        for aux in self.obsrun['frames']:
            filename = os.path.abspath(os.path.join(self.workenv['workdir'], aux))
            lista.append(DataFrame(filename=filename))

        self.obsresult = ObservationResult(mode='fiber_flat_image')
        self.obsresult.images=lista
        self.obsresult.instrument = 'MEGARA'

        self.requirement = [('ObservationResult',self.obsresult)]


    def __createTask(self, end, running, start):
        #TODO ver si es necesario dar el path completo. Puede estar en virtuenv o sepa dios donde...
        path = os.path.dirname(os.path.abspath(__file__))
        parameters = {'path':path,
                      'end':end,
                      'start':start,
                      'running':running}
        self.task = parameters


    def store(self):
        super(FiberFlatRecipe,self).store(self.result, self.__class__.__name__)

    def run(self):
        import time

        start = time.time()

        # TODO Desde aqui hasta "result = fits.HDUList([hdu, varhdu, num])" es el meodo run de cBias solo que con flow distinto
        cdata = []
        basicflow = SerialFlow(self.flow)

        try:
            for frame in self.obsresult.images:
                cdata.append(basicflow(frame.open()))

            data = c_median([d[0].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0][0].header)
        finally:
            for hdulist in cdata:
                hdulist.close()

        self.set_base_headers(hdu.header,  '0.14.dev0', 'FiberFlatRecipe', '0.1.0', 'FIBER_FLAT', 'MASTER_FIBER_FLAT', data[0].mean())

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        result = fits.HDUList([hdu, varhdu, num])

        _logger.info('prereduction ended')

        fiberflat_frame = result

        cstart = 2000
        step = 2

        # tracemap = self.trace(reduced[0].data, cstart, step)
        #
        # rss = apextract_tracemap(reduced[0].data, tracemap)
        #
        # rss[rss <= 0] = 1
        #
        # rss_norm = rss / rss.mean()
        #
        # _logger.info('fiber flat reduction ended')
        #
        # result = self.create_result(fiberflat_frame=reduced,
        #                             fiberflat_rss=fits.PrimaryHDU(rss_norm),
        #                             traces=tracemap)

        end = time.time()

        #TODO meter los otros dos resultados cuando no den error

        self.result.append(fiberflat_frame[0])

        self.__createTask(end, end - start, start)

        return result

    def trace(self, data, cstart, step):

        # fit_traces = domefun(data, cstart=2000, hs=20)

        cstart = cstart
        hs = 1
        step1 = step
        background1 = 10.0
        npred = 3
        maxdis1 = 2.0

        _logger.info('find peaks in column %i', cstart)

        central_peaks = init_traces(data, center=cstart, hs=hs,
                                background=background1, npred=npred)

        _logger.info(' %i peaks found', len(central_peaks))

        tracelist = []
        if data.dtype.byteorder != '=':
            _logger.debug('byteswapping image')
            image2 = data.byteswap().newbyteorder()
        else:
            image2 = data

        _logger.info('trace peaks')
        for trace in central_peaks.values():
            x, y, p = trace.start
            mm = trace(image2, x=x, y=y,
                         step=step1, hs=hs,
                         background=background1, maxdis=maxdis1
                         )

            pfit = numpy.polyfit(mm[:,0], mm[:,1], deg=5)

            tracelist.append({'fibid': trace.fibid, 'boxid': trace.boxid,
                              'start':0, 'stop':4095,
                              'fitparms': pfit.tolist()})

            #plt.title('cython version, fiber %i' % trace.fibid)
            #plt.plot(mm[:,0], mm[:,1], 'r*')
            #xpix = np.arange(0, 2000, 1)
            #p = numpy.poly1d(pfit)
            #plt.plot(xpix, p(xpix), 'b')
            #plt.show()

        return tracelist


# class RecipeInput(object):
#     def __init__(self, master_bias, obresult):
#         self.obresult = obresult
#         self.master_bias = master_bias
#
#         mi_clase = FiberFlatRecipe().run(obresult)
# rinput = RecipeInput(master_bias=DataFrame(filename='nombre'), obresult=None)
#
#
# ss = FiberFlatRecipe()
# obresult=None
# ss.run(rinput)



if __name__ == '__main__':

    workenv = {'datadir': '/home/pica/Documents/Megara/_data',
              'basedir': '/home/pica/Documents/Megara',
              'resultsdir': '/home/pica/Documents/Megara/_results',
              'workdir': '/home/pica/Documents/Megara/_work'}

    obsrun = [{'frames': ['r00002.fits'], 'instrument': 'MEGARA', 'pipeline': 'default', 'id': 1, 'mode': 'bias_image'}]

