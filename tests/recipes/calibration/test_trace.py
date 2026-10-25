"""The plots of the trace recipe are saved without pyplot"""

import matplotlib.pyplot as plt
import numpy

from megaradrp.recipes.calibration.trace import refine_boxes_from_array


def test_refine_boxes_from_array(tmp_path, monkeypatch):
    """The frontiers between the boxes are found 7 pixels after the expected ones"""
    monkeypatch.chdir(tmp_path)
    expected = [300, 600, 900, 1200, 1500]
    y = numpy.arange(1800.0)
    profile = 1000.0 + 50 * numpy.cos(y / 5)
    for frontier in expected:
        profile -= 800 * numpy.exp(-0.5 * ((y - frontier - 7) / 3) ** 2)
    arr = numpy.repeat(profile[:, numpy.newaxis], 20, axis=1)
    before = plt.get_fignums()
    refined, col = refine_boxes_from_array(arr, expected, 10, intermediate_results=True)
    assert refined == [frontier + 7 for frontier in expected]
    assert col == 10
    names = {path.name for path in tmp_path.iterdir()}
    assert names == {"offset_borders_cross.png", "frontiers_between_pseudoslits.png"}
    assert plt.get_fignums() == before
