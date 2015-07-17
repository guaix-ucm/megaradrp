__author__ = 'Pica4x6'

import os

from numina.core.oresult import ObservationResult
from numina.core.dataframe import DataFrame

from megaradrp.recipes.cBase import BaseRecipe


class DarkRecipe(BaseRecipe):
    '''Process BIAS images and create MASTER_BIAS.'''

    def __init__(self, obsrun , workenv):
        super(DarkRecipe,self).__init__(workenv)

        self.flow = []
        self.result = None
        self.obsrun = obsrun

        #TODO Poner esto mejor
        lista = []
        for aux in self.obsrun['frames']:
            filename = os.path.abspath(os.path.join(self.workenv['workdir'], aux))
            lista.append(DataFrame(filename=filename))
        self.obsresult = ObservationResult(mode='dark_image')
        self.obsresult.images=lista
        self.obsresult.instrument = 'MEGARA'

    @property
    def requirement(self):
        req = {'MasterBiasRequirement':''}
        return req

    def __createTask(self, end, running, start):
        #TODO ver si es necesario dar el path completo. Puede estar en virtuenv o sepa dios donde...
        path = os.path.dirname(os.path.abspath(__file__))
        parameters = {'path':path,
                      'end':end,
                      'start':start,
                      'running':running}
        self.task = parameters

    def store(self):
        super(DarkRecipe,self).store(self.result, self.__class__.__name__)


    def run(self, rinput=None):
        # TODO cDark. Metodo run()
        pass

