#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import json
import pathlib
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
    state.setdefault("flags", [])

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

    with NamedTemporaryFile(mode="w") as fd:
        json.dump(state, fd, cls=structured.ExtEncoder)
        fd.flush()  # To be sure the contents are in the file after dump
        my_open_file = tm.TraceMap._datatype_load(fd.name)

    assert my_open_file.instrument == state["instrument"]
    assert my_open_file.tags == state["tags"]
    assert my_open_file.uuid == state["uuid"]
    assert my_open_file.contents == state["contents"]
    assert my_open_file.DATATYPE == MegaraDataType.TRACE_MAP


def test_dump_tracemap(tracemap_data_state, benchmark=None):

    data, state = tracemap_data_state

    with NamedTemporaryFile() as my_file:
        # Using NamedTemporaryFile to generate a file name
        try:
            my_open_file = tm.TraceMap._datatype_dump(data, my_file.name)
            final = tm.TraceMap._datatype_load(my_open_file)
            traces = final.__getstate__()
        finally:
            pathlib.Path(my_open_file).unlink(missing_ok=True)

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


def test_geotrace_flags():
    """The flags are kept in the state, and some of them make the trace invalid"""
    mm = tm.GeometricTrace(fibid=10, boxid=1, start=6, stop=4090, fitparms=[100.0, 0.01])
    assert mm.flags == []
    assert mm.valid
    mm.add_flag("completed")
    mm.add_flag("completed")
    assert mm.flags == ["completed"]
    assert mm.valid
    mm.add_flag("crossing")
    assert not mm.valid
    state = mm.__getstate__()
    assert state["flags"] == ["completed", "crossing"]
    new = tm.GeometricTrace(**state)
    assert new.flags == ["completed", "crossing"]
    assert not new.valid


def test_setstate_geotrace_without_flags():
    """Traces saved before the flags have no flags"""
    new = tm.GeometricTrace.__new__(tm.GeometricTrace)
    new.__setstate__({"fibid": 1, "boxid": 0, "start": 6, "stop": 4090, "fitparms": [1.0]})
    assert new.flags == []
    assert new.valid


def test_geotrace_non_finite_is_not_valid():
    """A trace with coefficients that are not finite is not valid, also when loaded"""
    mm = tm.GeometricTrace(fibid=10, boxid=1, start=6, stop=4090, fitparms=[float("nan")] * 6)
    assert not mm.valid
    new = tm.GeometricTrace.__new__(tm.GeometricTrace)
    new.__setstate__({"fibid": 1, "boxid": 0, "start": 6, "stop": 4090, "fitparms": [100.0, float("inf")]})
    assert not new.valid
