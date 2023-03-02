from megaradrp.processing.trimover import OverscanCorrector
import numpy

def test_OverscanCorrector():

    detconf = {
        'trim1': [[0, 2056], [50, 4146]],
        'trim2': [[2156, 4212], [50, 4146]],
        'bng': [1, 1],
        'overscan1': [[0, 2056], [4149, 4196]],
        'overscan2': [[2156, 4212], [0, 50]],
        'prescan1': [[0, 2056], [0, 50]],
        'prescan2': [[2156, 4212], [4145, 4196]],
        'middle1': [[2056, 2106], [50, 4146]],
        'middle2': [[2106, 2156], [50, 4146]]
    }

    oc = OverscanCorrector(detconf)
    assert oc.trim1 == (slice(0, 2056, None),slice(50, 4146, None))
    assert oc.pcol1 == (slice(0, 2056, None),slice(0, 50, None))
    assert oc.ocol1 == (slice(0, 2056, None),slice(4149, 4196, None))
    assert oc.orow1 == (slice(2056, 2106, None),slice(50, 4146, None))

    assert oc.trim2 == (slice(2156, 4212, None),slice(50, 4146, None))
    assert oc.pcol2 == (slice(2156, 4212, None),slice(4145, 4196, None))
    assert oc.ocol2 == (slice(2156, 4212, None),slice(0, 50, None))
    assert oc.orow2 == (slice(2106, 2156, None),slice(50, 4146, None))


def test_amp_1():

    detconf = {
        'trim1': [[0, 2056], [50, 4146]],
        'trim2': [[2156, 4212], [50, 4146]],
        'bng': [1, 1],
        'overscan1': [[0, 2056], [4149, 4196]],
        'overscan2': [[2156, 4212], [0, 50]],
        'prescan1': [[0, 2056], [0, 50]],
        'prescan2': [[2156, 4212], [4145, 4196]],
        'middle1': [[2056, 2106], [50, 4146]],
        'middle2': [[2106, 2156], [50, 4146]]
    }

    oc = OverscanCorrector(detconf)

    data = 2000.0 + numpy.zeros(shape=(4212, 4196))

    # k,c,deg = spl1._eval_args

    fit1, spl1 = oc.eval_spline_amp1(data)
    data[oc.trim1] -= fit1[:, numpy.newaxis]

    assert fit1.shape == (2056,)

    fit2, spl2 = oc.eval_spline_amp2(data)
    data[oc.trim2] -= fit2[:, numpy.newaxis]

    assert fit2.shape == (2056,)


if __name__ == "__main__":
    test_OverscanCorrector()
