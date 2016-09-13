#
# Copyright 2015-2016 Universidad Complutense de Madrid
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


import megaradrp.products
from tempfile import NamedTemporaryFile
import pytest
import yaml


def create_test_tracemap():
    instrument = 'TEST1'
    tags = {}
    uuid = '123456789'
    data = megaradrp.products.TraceMap(instrument=instrument)
    data.tags = tags
    data.uuid = uuid

    state = dict(instrument=instrument,
                 tags=tags,
                 uuid=uuid,
                 tracelist=[]
                 )

    return data, state


def test_getstate_traceMap():

    data, state = create_test_tracemap()

    assert (data.__getstate__() == state)


def test_setstate_traceMap():

    data, state = create_test_tracemap()

    result = megaradrp.products.TraceMap(instrument='unknown')
    result.__setstate__(state)

    assert (state['instrument'] == result.instrument)
    assert (state['tags'] == result.tags)
    assert (state['uuid'] == result.uuid)
    assert (state['tracelist'] == result.tracelist)


@pytest.mark.xfail
def test_fail_traceMap():
    my_obj = megaradrp.products.TraceMap()
    my_obj._datatype_load('')


def test_load_traceMap():

    _data, state = create_test_tracemap()
    my_file = NamedTemporaryFile()

    with open(my_file.name, 'w') as fd:
        yaml.dump(state, fd)

    my_obj = megaradrp.products.TraceMap()
    my_open_file = my_obj._datatype_load(my_file.name)

    assert (my_open_file.instrument == state['instrument'])
    assert (my_open_file.tags == state['tags'])
    assert (my_open_file.uuid == state['uuid'])
    assert (my_open_file.tracelist == state['tracelist'])


def test_dump_traceMap(benchmark=None):

    class Aux(object):
        def __init__(self, destination):
            self.destination = destination

    data, state = create_test_tracemap()

    my_obj = megaradrp.products.TraceMap()
    my_file = NamedTemporaryFile()
    work_env = Aux(my_file.name)
    my_open_file = my_obj._datatype_dump(data, work_env)

    with open(my_open_file, 'r') as fd:
        traces = yaml.load(fd)

    assert (traces == state)


if __name__ == "__main__":
    test_load_traceMap()
    test_dump_traceMap()