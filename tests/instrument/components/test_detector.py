from megaradrp.testing.create_image import create_detector


def test_detector_shape():
    det = create_detector()
    img = det.readout()
    assert img.shape == (4212, 4196)
