
import pytest

import megaradrp.requirements as reqs
import numina.core
import numina.dal
from numina.exceptions import NoResultFound
from ..types import MasterFiberFlat
from ..products.tracemap import TraceMap


class IoInstance(object):
    _numina_desc_val = {}


class LocalDal(object):
    def __init__(self):
        self.filename = "some_file.txt"

    def search_product(self, name, tipo, obsres, options=None):
        if obsres.id == 1:
            return numina.dal.StoredProduct(1, self.filename, tags={})
        if obsres.id == 2:
            return numina.dal.StoredProduct(1, None, tags={})
        elif obsres.id == 0:
            raise NoResultFound('Nothing')
        else:
            sp = None
        return sp


@pytest.mark.parametrize('obsid', [0, 1, 2])
def test_extinction_none(obsid):

    name = "reference_extinction"
    reference_extinction = reqs.ReferenceExtinction()
    reference_extinction.dest = name

    dal = LocalDal()

    obs = numina.core.ObservationResult()
    obs.id = obsid
    obs.requirements = {name: None}
    value = reference_extinction.query(dal, obs)
    instance = IoInstance()

    reference_extinction.__set__(instance, value)
    cval = instance._numina_desc_val[name]
    assert cval is None


def _test_extinction_filename():

    name = "reference_extinction"
    reference_extinction = reqs.ReferenceExtinction()
    reference_extinction.dest = name

    dal = LocalDal()

    obs = numina.core.ObservationResult()

    import tempfile, os, numpy

    arr = numpy.array([1.0,2.0,3.0,4.0])
    fd, path = tempfile.mkstemp()
    try:
        #obs.requirements = {name: path}
        numpy.savetxt(path, arr)

        value = reference_extinction.query(dal, obs)
        instance = IoInstance()

        reference_extinction.__set__(instance, value)
        cval = instance._numina_desc_val[name]
        assert numpy.allclose(cval, arr)
        assert False
    finally:
        os.remove(path)


def test_tag1():
    import megaradrp.tests.simpleobj as simple

    fr = simple.create_simple_hdul()

    tipo = MasterFiberFlat()
    keys = ['instrument', 'uuid', 'observation_date']
    mu = tipo.extract_db_info(fr, keys)

    assert True


def test_tag2():
    import megaradrp.tests.simpleobj as simple

    tipo = TraceMap()
    fr = tipo
    keys = ['instrument', 'uuid', 'observation_date']
    mu = tipo.extract_db_info(fr, keys)

    assert True
