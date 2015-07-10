__author__ = 'Pica4x6'

import logging
import numpy

from astropy.io import fits

from numina.array.combine import median as c_median
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.flow import SerialFlow
from numina.flow.processing import BiasCorrector

from cBase import BaseRecipeAutoQC as MegaraBaseRecipe

from megaradrp.products import MasterFiberFlat
from megaradrp.products import TraceMap
from megaradrp.requirements import MasterBiasRequirement

from megaradrp.trace.traces import init_traces
from megaradrp.trace._traces import tracing  # @UnresolvedImport
from megaradrp.core import apextract2
from megaradrp.core import OverscanCorrector, TrimImage

_logger = logging.getLogger('numina.recipes.megara')


class FiberFlatRecipe(MegaraBaseRecipe):
    '''Process FIBER_FLAT images and create MASTER_FIBER_FLAT.'''

    # Requirements
    master_bias = MasterBiasRequirement()
    obresult = ObservationResultRequirement()
    # Products
    fiberflat_frame = Product(MasterFiberFlat)
    fiberflat_rss = Product(MasterFiberFlat)
    traces = Product(TraceMap)

    def __init__(self):
        super(FiberFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, rinput):
        return self.process_base1(rinput.obresult, rinput.master_bias)

    def process_common(self, obresult, master_bias):
        _logger.info('starting prereduction')

        o_c = OverscanCorrector()
        t_i = TrimImage()

        with master_bias.open() as hdul:
            mbias = hdul[0].data.copy()
            b_c = BiasCorrector(mbias)

        basicflow = SerialFlow([o_c, t_i, b_c])

        cdata = []

        try:
            for frame in obresult.frames:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                cdata.append(hdulist)

            _logger.info('stacking %d images using median', len(cdata))

            data = c_median([d[0].data for d in cdata], dtype='float32')
            template_header = cdata[0][0].header
            hdu = fits.PrimaryHDU(data[0], header=template_header)
        finally:
            for hdulist in cdata:
                hdulist.close()

        hdr = hdu.header
        hdr['IMGTYP'] = ('FIBER_FLAT', 'Image type')
        hdr['NUMTYP'] = ('MASTER_FIBER_FLAT', 'Data product type')
        hdr = self.set_base_headers(hdr)
        hdr['CCDMEAN'] = data[0].mean()

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        result = fits.HDUList([hdu, varhdu, num])

        _logger.info('prereduction ended')

        return result

    def process_base1(self, obresult, master_bias):
        _logger.info('starting fiber flat reduction')

        reduced = self.process_common(obresult, master_bias)

        cstart = 2000
        step = 2

        tracemap = self.trace(reduced[0].data, cstart, step)

        rss = apextract2(reduced[0].data, tracemap)

        rss[rss <= 0] = 1

        rss_norm = rss / rss.mean()

        _logger.info('fiber flat reduction ended')

        result = self.create_result(fiberflat_frame=reduced,
                                    fiberflat_rss=fits.PrimaryHDU(rss_norm),
                                    traces=tracemap)

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

            mm = tracing(image2, x=x, y=y,
                         p=p, step=step1, hs=hs,
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