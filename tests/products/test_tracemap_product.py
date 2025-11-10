#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import json
from tempfile import NamedTemporaryFile

import numpy
import numpy.polynomial.polynomial as nppol
import numina.types.structured as structured
import pytest

from megaradrp.datatype import MegaraDataType
import megaradrp.products.tracemap as tm


state1 = {"fibid": 100, "boxid": 12, "start": 1000, "stop": 2000}

state2 = {"fibid": 100, "boxid": 12, "start": 1000, "stop": 2000, "fitparms": [1, 2, 3]}

state3 = {"fibid": 100, "boxid": 12, "start": 1000, "stop": 2000, "fitparms": []}


@pytest.mark.parametrize("state", [state1, state2])
def test_geotrace(state):

    mm = tm.GeometricTrace(**state)
    assert isinstance(mm.polynomial, nppol.Polynomial)

    for key in state:
        assert state[key] == getattr(mm, key)

    if "fitparms" not in state:
        params = []
        coeff = [0.0]
    else:
        params = state["fitparms"]
        coeff = params

    assert mm.fitparms == params
    assert numpy.all(mm.polynomial.coef == coeff)


@pytest.mark.parametrize("state", [state1, state2])
def test_getstate_geotrace(state):

    mm = tm.GeometricTrace(**state)

    # Expected result
    if "fitparms" not in state:
        state["fitparms"] = []

    newstate = mm.__getstate__()

    for key in state:
        assert state[key] == getattr(mm, key)

    assert newstate == state


@pytest.mark.parametrize("state", [state2, state3])
def test_setstate_geotrace(state):
    new = tm.GeometricTrace.__new__(tm.GeometricTrace)
    new.__setstate__(state)

    for key in state:
        assert state[key] == getattr(new, key)

    if state["fitparms"]:
        coeff = state["fitparms"]
    else:
        coeff = [0.0]

    assert numpy.all(new.polynomial.coef == coeff)
    assert isinstance(new.polynomial, nppol.Polynomial)


def test_getstate_tracemap(tracemap_data_state):

    data, state = tracemap_data_state

    assert data.__getstate__() == state


def test_setstate_tracemap(tracemap_data_state):

    data, state = tracemap_data_state

    result = tm.TraceMap(instrument="unknown")
    result.__setstate__(state)

    assert state["instrument"] == result.instrument
    assert state["tags"] == result.tags
    assert state["uuid"] == result.uuid
    assert state["contents"] == result.contents


def test_fail_tracemap():
    with pytest.raises(FileNotFoundError):
        my_obj = tm.TraceMap()
        my_obj._datatype_load("")


def test_load_tracemap(tracemap_data_state):

    _data, state = tracemap_data_state
    my_file = NamedTemporaryFile()

    with open(my_file.name, "w") as fd:
        json.dump(state, fd, cls=structured.ExtEncoder)

    my_open_file = tm.TraceMap._datatype_load(my_file.name)

    assert my_open_file.instrument == state["instrument"]
    assert my_open_file.tags == state["tags"]
    assert my_open_file.uuid == state["uuid"]
    assert my_open_file.contents == state["contents"]
    assert my_open_file.DATATYPE == MegaraDataType.TRACE_MAP


def test_dump_tracemap(tracemap_data_state, benchmark=None):

    data, state = tracemap_data_state

    my_file = NamedTemporaryFile()
    my_open_file = tm.TraceMap._datatype_dump(data, my_file.name)

    final = tm.TraceMap._datatype_load(my_open_file)
    traces = final.__getstate__()

    assert traces == state


def test_query_fields():
    tracemap = tm.TraceMap()
    assert tracemap.query_expr.fields() == {"insmode", "vph"}
    assert tracemap.query_expr.tags() == {"insmode", "vph"}


def test_query_ds9(tmpdir, tracemap_data):
    p = tmpdir.join("ds9.reg")
    with p.open("w") as ds9reg:
        tracemap_data.to_ds9_reg(ds9reg)
    assert True
