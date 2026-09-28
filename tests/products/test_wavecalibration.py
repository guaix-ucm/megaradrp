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

import pytest
import numina.types.structured as structured

import megaradrp.products
import megaradrp.products.wavecalibration as wcal


def test_getstate_wavecalib(wavecalib_data_state):

    data, state = wavecalib_data_state

    assert data.__getstate__() == state


def test_setstate_wavecalib(wavecalib_data_state):

    data, state = wavecalib_data_state
    result = megaradrp.products.WavelengthCalibration(instrument="unknown")
    result.__setstate__(state)

    assert state["instrument"] == result.instrument
    assert state["tags"] == result.tags
    assert state["uuid"] == result.uuid
    for idx, cont in enumerate(state["contents"]):
        assert result.contents[idx].__getstate__() == cont


@pytest.mark.xfail
def test_fail_wcal():
    my_obj = megaradrp.products.WavelengthCalibration()
    my_obj._datatype_load("")


def test_load_wavecalib(wavecalib_data_state):

    _data, state = wavecalib_data_state

    with NamedTemporaryFile(mode="w") as fd:
        json.dump(state, fd, cls=structured.ExtEncoder)
        fd.flush()  # To be sure the contents are in the file after dump
        my_open_file = wcal.WavelengthCalibration._datatype_load(fd.name)

    assert my_open_file.instrument == state["instrument"]
    assert my_open_file.tags == state["tags"]
    assert my_open_file.uuid == state["uuid"]

    for idx, cont in enumerate(my_open_file.contents):
        assert cont.__getstate__() == state["contents"][idx]


def test_dump_wavecalib(wavecalib_data_state):

    data, state = wavecalib_data_state

    with NamedTemporaryFile() as my_file:
        # Using NamedTemporaryFile to generate a file name
        try:
            my_open_file = wcal.WavelengthCalibration._datatype_dump(data, my_file.name)
            final = wcal.WavelengthCalibration._datatype_load(my_open_file)
            traces = final.__getstate__()
        finally:
            pathlib.Path(my_open_file).unlink(missing_ok=True)

    assert traces == state


def test_query_fields():
    my_obj = megaradrp.products.WavelengthCalibration()
    assert my_obj.query_expr.fields() == {"insmode", "vph"}
    assert my_obj.query_expr.tags() == {"insmode", "vph"}
