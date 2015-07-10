__author__ = 'Pica4x6'

#
# Copyright 2008-2015 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

'''Basic tools and classes used to generate recipe modules.

A recipe is a class that complies with the *reduction recipe API*:

 * The class must derive from :class:`numina.core.BaseRecipe`.

'''

import abc
import traceback
import logging

from six import with_metaclass

from numina import __version__
from numina.core.recipeinout import ErrorRecipeResult, RecipeResult
from numina.core.recipeinout import RecipeRequirements
from numina.core.metarecipes import RecipeType,RecipeTypeAutoQC


_logger = logging.getLogger('numina')


class BaseRecipeMethods(object):
    '''Base class for all instrument recipes'''

    RecipeResult = RecipeResult
    RecipeRequirements = RecipeRequirements

    # Recipe own logger
    logger = _logger

    def __init__(self, *args, **kwds):
        super(BaseRecipeMethods, self).__init__()
        self.__author__ = 'Unknown'
        self.__version__ = '0.0.0'
        # These two are maintained
        # for the moment
        self.environ = {}
        self.runinfo = {}
        #
        self.instrument = None
        self.configure(**kwds)

    def configure(self, **kwds):
        if 'author' in kwds:
            self.__author__ = kwds['author']
        if 'version' in kwds:
            self.__version__ = kwds['version']
        if 'instrument' in kwds:
            self.instrument = kwds['instrument']
        if 'runinfo' in kwds:
            self.runinfo = kwds['runinfo']


    @classmethod
    def create_requirements(cls, *args, **kwds):
        '''
        Pass the result arguments to the RecipeRequirements constructor
        '''
        return cls.RecipeRequirements(*args, **kwds)

    @classmethod
    def create_result(cls, *args, **kwds):
        '''
        Pass the result arguments to the RecipeResult constructor
        '''
        return cls.RecipeResult(*args, **kwds)

    def run(self, recipe_input):
        return self.create_result()

    def __call__(self, recipe_input):
        '''
        Process the result of the observing block with the
        Recipe.

        :param recipe_input: the input appropriated for the Recipe
        :param type: RecipeRequirement
        :rtype: a RecipeResult object or an error

        '''

        try:
            result = self.run(recipe_input)
        except Exception as exc:
            _logger.error("During recipe execution %s", exc)
            return ErrorRecipeResult(
                exc.__class__.__name__,
                str(exc),
                traceback.format_exc()
                )
        return result

    def set_base_headers(self, hdr):
        '''Set metadata in FITS headers.'''
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        return hdr


class BaseRecipe(with_metaclass(abc.ABCMeta, BaseRecipeMethods)):
    '''Base class for all instrument recipes'''

    @abc.abstractmethod
    def run(self, recipe_input):
        return self.create_result()


class BaseRecipeAlt(with_metaclass(RecipeType, BaseRecipeMethods)):
    '''Base class for instrument recipes'''
    pass


class BaseRecipeAutoQC(with_metaclass(RecipeTypeAutoQC, BaseRecipeMethods)):
    '''Base class for instrument recipes'''
    pass
