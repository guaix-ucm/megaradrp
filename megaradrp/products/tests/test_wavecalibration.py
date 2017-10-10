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
    mm = SolutionArcCalibration(features, orig['coeff'], orig['residual_std'], cr_linear)
    return mm


orig = {
    'cr_linear': {
        'cdelt': 3.0, 'crmax': 4600, 'crmin': 2300, 'crpix': 1200, 'crval': 12},
    'features': [
        {'category': 'A', 'xpos': 100, 'reference': 3210, 'wavelength':1, 'line_ok': True, 'ypos': 0, 'peak': 0, 'funcost': 12.0, 'fwhm': 0,
         'lineid': 1},
        {'category': 'A', 'xpos': 150, 'reference': 3310, 'wavelength': 2, 'line_ok': True, 'ypos': 0, 'peak': 0, 'funcost': 12.0, 'fwhm': 0,
         'lineid': 2},
        {'category': 'C', 'xpos': 250, 'reference': 3410, 'wavelength': 3, 'line_ok': True, 'ypos': 0, 'peak': 0, 'funcost': 13.0, 'fwhm': 0,
         'lineid': 3}],
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
    data.error_fitting =  []
    data.missing_fibers = []
    data.total_fibers = 623
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
                 meta_info={},
                 type=data.name(),
                 contents=contents,
                 global_offset=[0.0],
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
def test_fail_traceMap():
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

    class Aux(object):
        def __init__(self, destination):
            self.destination = destination

    data, state = create_test_wavecalib()

    my_file = NamedTemporaryFile()
    work_env = Aux(my_file.name)
    my_open_file = wcal.WavelengthCalibration._datatype_dump(data, work_env)

    final = wcal.WavelengthCalibration._datatype_load(my_open_file)
    traces = final.__getstate__()

    assert (traces == state)


if __name__ == "__main__":
    test_dump_wavecalib()
    test_load_wavecalib()
