__author__ = 'Pica4x6'

import logging
_logger = logging.getLogger("numina")

class Megara(object):

    def __init__(self, directory, obsrun_list, requirement_file):
        self.workenv = directory
        self.obsrun = list(obsrun_list)
        self.requirement_file = requirement_file.copy()
        # TODO Hacerlo dinamicamente en funcion del nombre. Esta copiado de fichero megara.drp.yaml. Esta hecho por __init__
        self.recipes = {'bias_image':'megaradrp.recipes.calibration.cBias.BiasRecipe',
                        'dark_image': 'megaradrp.recipes.DarkRecipe',
                        'fiber_flat_image':'megaradrp.recipes.calibration.cFiberFlat.FiberFlatRecipe',
                        'mos_image':'megaradrp.recipes.scientific.FiberMOSRecipe2',
                        'flux_calibration':'megaradrp.recipes.calibration.PseudoFluxCalibrationRecipe',
                        'trace_map':'megaradrp.recipes.calibration.flat.TraceMapRecipe',
                        'fail':'numina.core.utils.AlwaysFailRecipe',
                        'success':'numina.core.utils.AlwaysSuccessRecipe'}

    @property
    def workenv(self):
        return self.__workenv

    @workenv.setter
    def workenv(self, directory):
        import os

        if not directory['workdir']:
            directory['workdir'] = os.path.abspath(os.path.join(directory['basedir'], '_work'))
            if not os.path.exists(directory['workdir']):
                os.makedirs(directory['workdir'])

        if not directory['resultsdir']:
            directory['resultsdir'] = os.path.abspath(os.path.join(directory['basedir'], '_results'))
            if not os.path.exists(directory['resultsdir']):
                os.makedirs(directory['resultsdir'])

        if not directory['datadir']:
            directory['datadir'] = os.path.abspath(os.path.join(directory['basedir'], '_data'))


        self.__workenv = {'basedir': directory['basedir'],
                        'workdir': directory['workdir'],
                        'resultsdir': directory['resultsdir'],
                        'datadir': directory['datadir']}

    def __deleteFrames(self, frames):
        import os
        import shutil
        shutil.rmtree(self.workenv['workdir'])
        os.makedirs(self.workenv['workdir'])

    def __copyFrames(self, frames):
        import shutil

        for frame in frames:
            ori = '%s/%s'%(self.workenv['datadir'],frame)
            dest = '%s/%s'%(self.workenv['workdir'],frame)
            shutil.copyfile(ori, dest)


    def run(self):
        import importlib

        for obs in self.obsrun:
            self.__copyFrames(obs['frames']+[self.requirement_file['requirements']['master_bias']])
            class_path = '.'.join(self.recipes[obs['mode']].split('.')[:-1])
            class_name = self.recipes[obs['mode']].split('.')[-1]
            module = importlib.import_module(class_path)
            class_ = getattr(module, class_name)
            #TODO hay que ver como afecta el fichero control.yaml=requirement_file a cada una de las observaciones
            instance = class_(obs, self.workenv, self.requirement_file)
            instance.run()
            instance.store()
            self.__deleteFrames(obs['frames'])


if __name__ == '__main__':
    test = Megara({'basedir':'/home/pica/Documents/Megara', 'workdir':'','resultsdir':'','datadir':''},{})
    print test.workenv
    print test.obsrun
