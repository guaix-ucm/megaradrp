#
# Copyright 2015-2017 Universidad Complutense de Madrid
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


from tempfile import NamedTemporaryFile
import json

import pytest
import numpy
import numpy.polynomial.polynomial as nppol

import megaradrp.products
import megaradrp.products.tracemap as tm


def create_test_tracemap():
    instrument = 'TEST1'
    tags = {}
    uuid = '123456789'
    data = megaradrp.products.TraceMap(instrument=instrument)
    data.tags = tags
    data.uuid = uuid
    data.total_fibers = 623
    state = dict(instrument=instrument,
                 tags=tags,
                 uuid=uuid,
                 error_fitting=[],
                 missing_fibers=[],
                 total_fibers=623,
                 meta_info={},
                 contents=[],
                 boxes_positions=[],
                 type=data.name(),
                 ref_column=2000,
                 global_offset=[0.0]
                 )

    return data, state


state1 = {
    'fibid': 100,
    'boxid': 12,
    'start': 1000,
    'stop': 2000
}

state2 = {
    'fibid': 100,
    'boxid': 12,
    'start': 1000,
    'stop': 2000,
    'fitparms': [1, 2, 3]
}

state3 = {
    'fibid': 100,
    'boxid': 12,
    'start': 1000,
    'stop': 2000,
    'fitparms': []
}

@pytest.mark.parametrize("state", [state1, state2])
def test_geotrace(state):

    mm = tm.GeometricTrace(**state)
    assert isinstance(mm.polynomial, nppol.Polynomial)

    for key in state:
        assert state[key] == getattr(mm, key)

    if 'fitparms' not in state:
        params = []
        coeff = [0.0]
    else:
        params = state['fitparms']
        coeff = params

    assert mm.fitparms == params
    assert numpy.all(mm.polynomial.coef == coeff)


@pytest.mark.parametrize("state",[state1, state2])
def test_getstate_geotrace(state):

    mm = tm.GeometricTrace(**state)

    # Expected result
    if 'fitparms' not in state:
        state['fitparms'] = []

    newstate = mm.__getstate__()

    for key in state:
        assert state[key] == getattr(mm, key)

    assert newstate == state


@pytest.mark.parametrize("state",[state2, state3])
def test_setstate_geotrace(state):
    new = tm.GeometricTrace.__new__(tm.GeometricTrace)
    new.__setstate__(state)

    for key in state:
        assert state[key] == getattr(new, key)

    if state['fitparms']:
        coeff = state['fitparms']
    else:
        coeff = [0.0]

    assert numpy.all(new.polynomial.coef == coeff)
    assert isinstance(new.polynomial, nppol.Polynomial)


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
    assert (state['contents'] == result.contents)


@pytest.mark.xfail
def test_fail_traceMap():
    my_obj = megaradrp.products.TraceMap()
    my_obj._datatype_load('')


def test_load_traceMap():

    _data, state = create_test_tracemap()
    my_file = NamedTemporaryFile()

    with open(my_file.name, 'w') as fd:
        json.dump(state, fd)

    my_obj = megaradrp.products.TraceMap()
    my_open_file = my_obj._datatype_load(my_file.name)

    assert (my_open_file.instrument == state['instrument'])
    assert (my_open_file.tags == state['tags'])
    assert (my_open_file.uuid == state['uuid'])
    assert (my_open_file.contents == state['contents'])


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
        traces = json.load(fd)

    assert (traces == state)


if __name__ == "__main__":
    test_load_traceMap()
    test_dump_traceMap()
