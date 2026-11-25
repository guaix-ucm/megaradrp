"""Quality control of the model map"""

import numpy
from numina.types.qc import QC

from megaradrp.processing.modelmap import check_model_map
from megaradrp.products.modelmap import GeometricModel, ModelMap

NBOX, PERBOX, PITCH, GAP = 3, 10, 6.0, 14.0
NCOL = 4096


def fiber_centre(index):
    box, pos = divmod(index, PERBOX)
    return 100.0 + box * (PERBOX * PITCH + GAP) + pos * PITCH


def make_model_map(centre=None, width=None):
    """3 blocks of 10 parallel fibers, with widths that grow towards the borders

    `centre` and `width` are optional functions (fibid, x) -> correction.
    """
    model_map = ModelMap(instrument="MEGARA")
    model_map.total_fibers = NBOX * PERBOX
    for index in range(NBOX * PERBOX):
        fibid = index + 1
        y0 = fiber_centre(index)

        def mean(x, y0=y0, fibid=fibid):
            x = numpy.asarray(x, dtype=float)
            out = y0 + 1e-3 * x + 2e-7 * (x - 2000) ** 2
            return out + (centre(fibid, x) if centre else 0.0)

        def stddev(x, y0=y0, fibid=fibid):
            x = numpy.asarray(x, dtype=float)
            out = 1.5 + 2e-8 * (x - 2000) ** 2 + 1e-4 * abs(y0 - 200)
            return out + (width(fibid, x) if width else 0.0)

        model = {"model_name": "gaussbox", "params": {"mean": mean, "stddev": stddev}}
        model_map.contents.append(GeometricModel(fibid, index // PERBOX + 1, 1, NCOL, model))
    return model_map


def flat_rss(model_map, change=None):
    """RSS of a flat: a lamp spectrum with absorption bands, a throughput per fiber"""
    x = numpy.arange(NCOL)
    lamp = 1000.0 * (1 + 0.3 * numpy.sin(x / 800.0))
    lamp[1500:1520] *= 0.7
    rng = numpy.random.default_rng(3)
    throughput = rng.uniform(0.5, 1.2, model_map.total_fibers)
    rss = throughput[:, numpy.newaxis] * lamp * (1 + 0.05 * x / NCOL)
    if change:
        change(rss)
    return rss


def flags(model_map):
    return {m.fibid: m.flags for m in model_map.contents if m.flags}


def test_good():
    model_map = make_model_map()
    assert check_model_map(model_map, flat_rss(model_map)) == {}
    assert model_map.quality_control == QC.GOOD
    assert flags(model_map) == {}


def test_not_parallel():
    """A fiber not parallel to its neighbours; the neighbours are not flagged"""
    model_map = make_model_map(centre=lambda fibid, x: 3e-4 * x if fibid == 15 else 0.0)
    found = check_model_map(model_map)
    assert found == {"not_parallel": [15]}
    assert flags(model_map) == {15: ["not_parallel"]}
    assert model_map.quality_control == QC.PARTIAL


def test_not_parallel_border():
    """A fiber at the border of a block is compared with its only neighbour"""
    model_map = make_model_map(centre=lambda fibid, x: 3e-4 * x if fibid == 11 else 0.0)
    assert check_model_map(model_map) == {"not_parallel": [11]}


def test_crossing():
    """Two traces that cross; with a large deviation they are also not parallel"""
    model_map = make_model_map(centre=lambda fibid, x: 2.4e-3 * x if fibid == 5 else 0.0)
    found = check_model_map(model_map)
    assert found["crossing"] == [5, 6]
    assert 5 in found["not_parallel"]


def test_width():
    model_map = make_model_map(width=lambda fibid, x: 0.6 * (x > 2500) if fibid == 25 else 0.0)
    assert check_model_map(model_map) == {"width_not_smooth": [25]}


def test_amplitude():
    """A local change of the flux of one fiber, not of its neighbours; the lamp bands are common"""

    def dip(rss):
        rss[7, 3000:3150] *= 0.8

    model_map = make_model_map()
    assert check_model_map(model_map, flat_rss(model_map, dip)) == {"amplitude_not_smooth": [8]}


def test_bad():
    """More than 5 % of the fibers flagged"""
    model_map = make_model_map(centre=lambda fibid, x: 3e-4 * x if fibid in (4, 15, 26) else 0.0)
    assert check_model_map(model_map) == {"not_parallel": [4, 15, 26]}
    assert model_map.quality_control == QC.BAD


def test_invalid_models_are_skipped():
    model_map = make_model_map()
    model_map.contents[3].model = {}
    assert check_model_map(model_map) == {}
    assert model_map.quality_control == QC.GOOD


def test_no_models():
    model_map = ModelMap(instrument="MEGARA")
    assert check_model_map(model_map) == {}
    assert model_map.quality_control == QC.BAD
