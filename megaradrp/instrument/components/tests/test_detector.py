
from ..detector import MegaraDetector, ReadParams


def create_detector():
    DSHAPE = (2056 * 2, 2048 * 2)
    PSCAN = 50
    OSCAN = 50
    qe = 1.0
    dcurrent = 3.0 / 3600

    readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000.0)
    readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005.0)

    detector = MegaraDetector(
        'Detector',
        DSHAPE, OSCAN, PSCAN, qe=qe, dark=dcurrent,
        readpars1=readpars1, readpars2=readpars2, bins='11'
    )
    return detector


def test_detector_shape():
    det = create_detector()
    img = det.readout()
    assert img.shape == (4212, 4196)