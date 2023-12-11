#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import megaradrp.products
from tempfile import NamedTemporaryFile
import pytest
import json

import numina.types.qc
import numina.types.structured as structured
from numina.array.wavecalib.arccalibration import SolutionArcCalibration, WavecalFeature, CrLinear

import megaradrp.products.wavecalibration as wcal


# FIXME: copied from numina, this should be a fixture
def create_solution(orig):

    # Create Features
    features = []
    for feature in orig['features']:
        features.append(
            WavecalFeature(**feature)
        )

    cr_linear = CrLinear(**orig['cr_linear'])
    mm = SolutionArcCalibration(
        features, orig['coeff'], orig['residual_std'], cr_linear)
    return mm


orig = {
    'cr_linear': {
        'cdelt': 3.0, 'crmax': 4600, 'crmin': 2300, 'crpix': 1200, 'crval': 12},
    'features': [
        {'category': 'A', 'xpos': 100, 'reference': 3210, 'wavelength': 1, 'line_ok': True, 'ypos': 0,
         'peak': 0, 'funcost': 12.0, 'fwhm': 0, 'lineid': 1},
        {'category': 'A', 'xpos': 150, 'reference': 3310, 'wavelength': 2, 'line_ok': True, 'ypos': 0,
         'peak': 0, 'funcost': 12.0, 'fwhm': 0, 'lineid': 2},
        {'category': 'C', 'xpos': 250, 'reference': 3410, 'wavelength': 3, 'line_ok': True, 'ypos': 0,
         'peak': 0, 'funcost': 13.0, 'fwhm': 0, 'lineid': 3}],
    'coeff': [1.0, 0.1],
    'residual_std': 1.0,

}


def create_test_wavecalib():
    instrument = 'TEST1'
    tags = {'insmode': 'LCB', 'vph': 'LR-I'}
    uuid = '123456789'
    data = wcal.WavelengthCalibration(instrument=instrument)
    data.tags = tags
    data.uuid = uuid
    data.error_fitting = []
    data.missing_fibers = []
    data.total_fibers = 623
    meta_info = wcal.WavelengthCalibration.create_meta_info()
    meta_info['instrument_name'] = instrument
    meta_info['creation_date'] = data.meta_info['creation_date']
    contents = []

    for key in range(10):
        fibid = key + 1
        solution = create_solution(orig)
        fibersolution = wcal.FiberSolutionArcCalibration(fibid, solution)
        data.contents.append(fibersolution)
        contents.append(data.contents[-1].__getstate__())

    for key in range(100, 101):
        fibid = key + 1
        solution = create_solution(orig)
        fibersolution = wcal.FiberSolutionArcCalibration(fibid, solution)
        data.contents.append(fibersolution)
        contents.append(data.contents[-1].__getstate__())

    state = dict(instrument=instrument,
                 tags=tags,
                 uuid=uuid,
                 total_fibers=623,
                 missing_fibers=[],
                 error_fitting=[],
                 meta_info=meta_info,
                 type=data.name(),
                 contents=contents,
                 global_offset=[0.0],
                 type_fqn='megaradrp.products.wavecalibration.WavelengthCalibration',
                 quality_control=numina.types.qc.QC.UNKNOWN
                 )
    return data, state


def test_getstate_wavecalib():

    data, state = create_test_wavecalib()

    assert (data.__getstate__() == state)


def test_setstate_wavecalib():

    data, state = create_test_wavecalib()
    result = megaradrp.products.WavelengthCalibration(instrument='unknown')
    result.__setstate__(state)

    assert (state['instrument'] == result.instrument)
    assert (state['tags'] == result.tags)
    assert (state['uuid'] == result.uuid)
    for idx, cont in enumerate(state['contents']):
        assert (result.contents[idx].__getstate__() == cont)


@pytest.mark.xfail
def test_fail_wcal():
    my_obj = megaradrp.products.WavelengthCalibration()
    my_obj._datatype_load('')


def test_load_wavecalib():

    _data, state = create_test_wavecalib()
    my_file = NamedTemporaryFile()

    with open(my_file.name, 'w') as fd:
        json.dump(state, fd, cls=structured.ExtEncoder)

    my_open_file = wcal.WavelengthCalibration._datatype_load(my_file.name)

    assert (my_open_file.instrument == state['instrument'])
    assert (my_open_file.tags == state['tags'])
    assert (my_open_file.uuid == state['uuid'])

    for idx, cont in enumerate(my_open_file.contents):
        assert (cont.__getstate__() == state['contents'][idx])


def test_dump_wavecalib():

    data, state = create_test_wavecalib()

    my_file = NamedTemporaryFile()
    my_open_file = wcal.WavelengthCalibration._datatype_dump(
        data, my_file.name)

    final = wcal.WavelengthCalibration._datatype_load(my_open_file)
    traces = final.__getstate__()

    assert (traces == state)


def test_query_fields():
    my_obj = megaradrp.products.WavelengthCalibration()
    assert my_obj.query_expr.fields() == {'insmode', 'vph'}
    assert my_obj.query_expr.tags() == {'insmode', 'vph'}


if __name__ == "__main__":
    test_dump_wavecalib()
    test_load_wavecalib()
