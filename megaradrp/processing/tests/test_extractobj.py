
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
