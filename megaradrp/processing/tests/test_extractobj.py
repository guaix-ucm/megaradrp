
import numpy

import megaradrp.processing.extractobj as eobj


def test_cover_1():

    arr = [1, 1, 1, 1, 1, 1]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 7)

    arr = [0, 0, 0, 0, 0, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 0)


def test_cover_2():

    arr = [0, 1, 1, 1, 0, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 4, None)

    arr = [0, 1, 1, 1, 1, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 5)

    arr = [0, 1, 1, 1, 1, 1]
    res = eobj.coverage_det(arr)
    assert res == slice(1, 7)

    arr = [1, 1, 1, 1, 1, 0]
    res = eobj.coverage_det(arr)
    assert res == slice(0, 5)


def test_coverage_issue247():

    # Example from test case
    fmap = numpy.zeros((4300,), dtype='int32')
    fmap[226:4136] = 37
    fmap[224] = 1
    fmap[225] = 10
    fmap[4136] = 32

    # Test total coverage
    arr1 = fmap == fmap.max()
    res1 = eobj.coverage_det(arr1.astype('int'))
    assert res1 == slice(226, 4136, None)

    # Test max coverage
    arr2 = fmap > 0
    res2 = eobj.coverage_det(arr2.astype('int'))
    assert res2 == slice(224, 4137, None)
