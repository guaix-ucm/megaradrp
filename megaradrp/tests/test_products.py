#
# Copyright 2015 Universidad Complutense de Madrid
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

from megaradrp.products import TraceMap
from tempfile import NamedTemporaryFile
import pytest
import yaml

@pytest.mark.xfail
def test_fail_traceMap(benchmark=None):
    my_obj = TraceMap()
    my_obj._datatype_load('')

def test_load_traceMap(benchmark=None):
    data = dict(A = 'a',
                B = dict(C = 'c',
                         D = 'd',
                         E = 'e',)
                )

    my_obj = TraceMap()
    my_file = NamedTemporaryFile()
    with open(my_file.name, 'w') as fd:
        yaml.dump(data, fd)
    my_open_file = my_obj._datatype_load(my_file.name)
    assert (my_open_file == {'A': 'a', 'B': {'C': 'c', 'E': 'e', 'D': 'd'}})


def test_dump_traceMap(benchmark=None):
    class aux(object):
        def __init__(self, destination):
            self.destination = destination

    data = dict(A = 'a',
                B = dict(C = 'c',
                         D = 'd',
                         E = 'e',)
                )

    my_obj = TraceMap()
    my_file = NamedTemporaryFile()
    work_env = aux(my_file.name)
    my_open_file = my_obj._datatype_dump(data, work_env)
    with open(my_open_file, 'r') as fd:
        traces = yaml.load(fd)
    assert (traces == {'A': 'a', 'B': {'C': 'c', 'E': 'e', 'D': 'd'}})

if __name__ == "__main__":
    test_load_traceMap()
    test_dump_traceMap()