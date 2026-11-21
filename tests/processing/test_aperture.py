import numpy

from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.products.modelmap import ModelMap
from megaradrp.testing.create_traceoffset import create_offset_tracemap


def test_aperture_extractor_offset_copy():
    """The offset is applied to a copy: the calibration keeps its global offset"""
    tracemap = create_offset_tracemap()
    tracemap.global_offset = numpy.polynomial.Polynomial([0.5])
    extractor = ApertureExtractor(tracemap, offset=[1.5, 1e-3])
    assert numpy.allclose(tracemap.global_offset.coef, [0.5])
    assert numpy.allclose(extractor.trace_repr.global_offset.coef, [2.0, 1e-3])


def test_aperture_extractor_offset_modelmap_matrices():
    """The extraction matrices of a ModelMap, computed with the old offset, are not reused"""
    modelmap = ModelMap(instrument="MEGARA")
    modelmap._wcols = "matrices of the old offset"
    extractor = ApertureExtractor(modelmap, offset=[1.0])
    assert extractor.trace_repr._wcols is None
    assert modelmap._wcols == "matrices of the old offset"
