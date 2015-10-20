#
# Copyright 2011-2014 Universidad Complutense de Madrid
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


from megaradrp.taggers import get_tags_from_full_ob,tagger_empty,tagger_vph
from numina.core.oresult import obsres_from_dict
import yaml
import os, inspect
import pytest

def test_get_tags_from_full_ob(benchmark):
    loaded_obs = {}
    with open(os.path.dirname(inspect.getfile(inspect.currentframe()))+'/obsrun_image.yaml') as fd:
        for doc in yaml.load_all(fd):
            loaded_obs[doc['id']] = doc

    obsres = obsres_from_dict((loaded_obs[doc['id']]))
    obsres.frames[0].filename = os.path.dirname(inspect.getfile(inspect.currentframe()))+'/test.fits'
    obj = benchmark(get_tags_from_full_ob,obsres)
    assert obj == {}



def test_tagger_empty(benchmark):
    loaded_obs = {}
    with open(os.path.dirname(inspect.getfile(inspect.currentframe()))+'/obsrun_image.yaml') as fd:
        for doc in yaml.load_all(fd):
            loaded_obs[doc['id']] = doc

    obsres = obsres_from_dict((loaded_obs[doc['id']]))
    obsres.frames[0].filename = os.path.dirname(inspect.getfile(inspect.currentframe()))+'/test.fits'
    obj = benchmark(tagger_empty,obsres)
    assert obj == {}



def test_tagger_vph(benchmark):
    loaded_obs = {}
    with open(os.path.dirname(inspect.getfile(inspect.currentframe()))+'/obsrun_image.yaml') as fd:
        for doc in yaml.load_all(fd):
            loaded_obs[doc['id']] = doc

    obsres = obsres_from_dict((loaded_obs[doc['id']]))
    obsres.frames[0].filename = os.path.dirname(inspect.getfile(inspect.currentframe()))+'/test.fits'
    obj = benchmark(tagger_vph,obsres)
    assert obj == {'vph': 'VPH405_LR'}


def test_get_tags_from_full_ob_no_images(benchmark):
    loaded_obs = {}
    with open(os.path.dirname(inspect.getfile(inspect.currentframe()))+'/obsrun_no_images.yaml') as fd:
        for doc in yaml.load_all(fd):
            loaded_obs[doc['id']] = doc

    obsres = obsres_from_dict((loaded_obs[doc['id']]))
    obj = benchmark(tagger_vph,obsres)
    assert obj == {}


if __name__ == "__main__":
    test_get_tags_from_full_ob()