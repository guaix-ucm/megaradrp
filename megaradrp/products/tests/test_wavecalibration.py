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

from numina.array.wavecalib.arccalibration import SolutionArcCalibration, WavecalFeature, CrLinear


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


orig = {'cr_linear': {'cdelt': 3.0, 'crmax': 4600, 'crmin': 2300, 'crpix': 1200,
                          'crval': 12},
            'features': [
        {'category': 'A', 'xpos': 100, 'wv': 3210, 'line_ok': True, 'ypos': 0, 'flux': 0, 'funcost': 12.0, 'fwhm': 0,
         'lineid': 1},
        {'category': 'A', 'xpos': 150, 'wv': 3310, 'line_ok': True, 'ypos': 0, 'flux': 0, 'funcost': 12.0, 'fwhm': 0,
         'lineid': 2},
        {'category': 'C', 'xpos': 250, 'wv': 3410, 'line_ok': True, 'ypos': 0, 'flux': 0, 'funcost': 13.0, 'fwhm': 0,
         'lineid': 3}], 'wavelength': [11.0, 16.0, 26.0], 'coeff': [1.0, 0.1], 'residual_std': 1.0}


def create_test_wavecalib():
    instrument = 'TEST1'
    tags = {}
    uuid = '123456789'
    data = megaradrp.products.WavelengthCalibration(instrument=instrument)
    data.tags = tags
    data.uuid = uuid

    contents = {}
    for key in range(10):
        data.contents[key] = create_solution(orig)
        contents[key] = data.contents[key].__getstate__()

    for key in range(100, 101):
        data.contents[key] = create_solution(orig)
        contents[key] = data.contents[key].__getstate__()

    state = dict(instrument=instrument,
                 tags=tags,
                 uuid=uuid,
                 contents=contents
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
    for key in state['contents']:
        assert (result.contents[key].__getstate__() == state['contents'][key])


@pytest.mark.xfail
def test_fail_traceMap():
    my_obj = megaradrp.products.WavelengthCalibration()
    my_obj._datatype_load('')


def test_load_wavecalib():

    _data, state = create_test_wavecalib()
    my_file = NamedTemporaryFile()

    with open(my_file.name, 'w') as fd:
        yaml.dump(state, fd)

    my_obj = megaradrp.products.WavelengthCalibration()
    my_open_file = my_obj._datatype_load(my_file.name)

    assert (my_open_file.instrument == state['instrument'])
    assert (my_open_file.tags == state['tags'])
    assert (my_open_file.uuid == state['uuid'])
    for key in my_open_file.contents:
        assert (my_open_file.contents[key].__getstate__() == state['contents'][key])


def test_dump_wavecalib():

    class Aux(object):
        def __init__(self, destination):
            self.destination = destination

    data, state = create_test_wavecalib()

    my_obj = megaradrp.products.WavelengthCalibration()
    my_file = NamedTemporaryFile()
    work_env = Aux(my_file.name)
    my_open_file = my_obj._datatype_dump(data, work_env)

    with open(my_open_file, 'r') as fd:
        traces = yaml.load(fd)

    assert (traces == state)


if __name__ == "__main__":
    test_dump_wavecalib()
    test_load_wavecalib()
