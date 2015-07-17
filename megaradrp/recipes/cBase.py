# __author__ = 'Pica4x6'
#
# #
# # Copyright 2008-2015 Universidad Complutense de Madrid
# #
# # This file is part of Numina
# #
# # Numina is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # Numina is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
# #
# # You should have received a copy of the GNU General Public License
# # along with Numina.  If not, see <http://www.gnu.org/licenses/>.

import abc
from numina.core.oresult import ObservationResult


class Requirement(list):
    '''
        Class Requirement is a double linked circular list
    '''
    def __init__(self, sequence=[]):
        super(Requirement, self).__init__(sequence)
        self.position = 0

    def current(self):
        return self[self.position]

    def next(self, n=1):
        self.position = (self.position + n) % len(self)
        return self[self.position]

    def prev(self, n=1):
        return self.next(-n)


class BaseRecipe(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, workenv, obsrun, requirement_file):
        self.workenv = workenv
        self.result = None
        self.obsrun = obsrun
        self.requirement_file = requirement_file


    @property
    def requirement(self):
        return self.__requirement

    @requirement.setter
    def requirement(self, lista_requirement):
        self.__requirement = Requirement(lista_requirement)

    @property
    def result(self):
        return self.__result

    @result.setter
    def result(self, result):
        self.__result = result

    @property
    def task(self):
        return self.__task

    @task.setter
    def task(self, parameters):
        # self.__task = _task.copy()

        print parameters

        self.__task = {'observation':{'instrument':self.obsrun['instrument'],
                                    'mode':self.obsrun['mode'],
                                    'observing_result':self.obsrun['id']},
                     'result':'result.yaml',
                     'runinfo':{'data_dir':self.workenv['datadir'],
                                'pipeline':self.obsrun['pipeline'],
                                'recipe':self.__class__.__name__,
                                'recipe_fulpathfile':parameters['path'],
                                'recipe_version':'0.1.0',
                                'results_dir':self.workenv['resultsdir'],
                                'runner':'numina',
                                'runner_version':'0.14.dev0',
                                'time_end':parameters['end'],
                                'time_running':parameters['running'],
                                'time_start':parameters['start'],
                                'work_dir':self.workenv['workdir']}
                     }

    @abc.abstractmethod
    def run(self):
        pass

    def getGTC(self):
        pass

    def store(self, hdulist, file_yaml):
        self.__storeHDU(hdulist, file_yaml)
        self.__storeYaml(file_yaml)

    def __storeHDU(self, hdulist, filename):
        for hdu in hdulist:
            hdu.writeto('%s/%s.fits'%(self.workenv['resultsdir'],filename), clobber=True)

    def __storeYaml(self, filename):
        import yaml

        with open('%s/Task_%s.yaml'%(self.workenv['resultsdir'],filename), 'w') as yaml_file:
            yaml_file.write( yaml.dump(self.task, default_flow_style=False))

    def set_base_headers(self, hdr, numxver, name, numrver, img_type, num_typ, ccd_mean ):
        '''Set metadata in FITS headers.'''
        # TODO Esto no me gusta asi. Tiene que hacer algo en algun lado para poner las cabeceras por referencia...
        hdr['NUMXVER'] = (numxver, 'Numina package version')
        hdr['NUMRNAM'] = (name, 'Numina recipe name')
        hdr['NUMRVER'] = (numrver, 'Numina recipe version')
        hdr['IMGTYP'] = (img_type, 'Image type')
        hdr['NUMTYP'] = (num_typ, 'Data product type')
        hdr['CCDMEAN'] = ccd_mean

        return hdr

