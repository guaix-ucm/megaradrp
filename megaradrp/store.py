#
# Copyright 2011-2015 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

import logging

import yaml

from numina.user.dump import dump
from numina.user.load import load

from .products import TraceMap

_logger = logging.getLogger('megaradrp')


_logger.debug('register dump functions')


@dump.register(TraceMap)
def _d(tag, obj, where):

    filename = where.destination + '.yaml'

    with open(filename, 'w') as fd:
        yaml.dump(obj, fd)

    return filename

_logger.debug('register load functions')


@load.register(TraceMap)
def _l(tag, obj):

    with open(obj, 'r') as fd:
        traces = yaml.load(fd)

    return traces
