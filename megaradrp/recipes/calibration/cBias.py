__author__ = 'Pica4x6'

import os

from astropy.io import fits

from megaradrp.processing import OverscanCorrector, TrimImage
from numina.core.oresult import ObservationResult
from numina.core.dataframe import DataFrame

from numina.core import RecipeError
from numina.array.combine import median as c_median
from numina.flow import SerialFlow

from megaradrp.recipes.cBase import BaseRecipe


class BiasRecipe(BaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    def __init__(self, obsrun , workenv, requirement_file):
        super(BiasRecipe,self).__init__(workenv, obsrun, requirement_file)

        self.flow = [OverscanCorrector(), TrimImage()]
        self.result = []

        #TODO Poner esto mejor y hacerlo con los requirement
        lista = []
        for aux in self.obsrun['frames']:
            filename = os.path.abspath(os.path.join(self.workenv['workdir'], aux))
            lista.append(DataFrame(filename=filename))

        self.obsresult = ObservationResult(mode='bias_image')
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
        super(BiasRecipe,self).store(self.result, self.__class__.__name__)


    def run (self):
        import time

        start = time.time()

        if not self.obsresult.images:
            raise RecipeError('Frame list is empty')

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

        self.set_base_headers(hdu.header,  '0.14.dev0', 'BiasRecipe', '0.1.0', 'BIAS', 'MASTER_BIAS', data[0].mean())

        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        hdulist = fits.HDUList([hdu, varhdu, num])

        self.result.append(hdulist[0])

        end = time.time()

        self.__createTask(end, end - start, start)


        # TODO ver que acepta GTC y en funcion de eso crear una nueva clase con type() o devolver json

        return hdulist



if __name__ == '__main__':

    workenv = {'datadir': '/home/pica/Documents/Megara/_data',
              'basedir': '/home/pica/Documents/Megara',
              'resultsdir': '/home/pica/Documents/Megara/_results',
              'workdir': '/home/pica/Documents/Megara/_work'}

    obsrun = [{'frames': ['r00002.fits'], 'instrument': 'MEGARA', 'pipeline': 'default', 'id': 1, 'mode': 'bias_image'}]

    for obs in obsrun:
        instance = BiasRecipe(obs, workenv)
        print instance.requirement
        # print "TUPLA: ", [item for item in instance.requirement if 'ObservationResult' in item]
        # instance.run()
        # instance.store()

        hdulist = fits.open('/home/pica/Documents/Megara/_results/biasframe.fits')
        print hdulist.info()
        prihdr = hdulist[0].header
        print prihdr.keys()