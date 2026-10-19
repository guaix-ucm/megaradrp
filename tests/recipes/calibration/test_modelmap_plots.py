"""The plots of the parameters of the fibers of the model map"""

import numpy as np
import pytest

from megaradrp.recipes.calibration.modelmap import save_parameter_plots


@pytest.mark.parametrize("processes", [1, 2])
def test_save_parameter_plots(tmp_path, monkeypatch, processes):
    monkeypatch.chdir(tmp_path)
    g_col = np.arange(100, 4100, 100)
    plots = [
        (f"fib_{fibid:03d}_{name}.png", f"{name} fib{fibid:03d}", 1.5 + 0.0 * g_col, 1.5 + 0.0 * g_col)
        for fibid in range(1, 6)
        for name in ["mean", "stddev"]
    ]
    save_parameter_plots(g_col, plots, processes=processes)
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == sorted(plot[0] for plot in plots)
    assert all((tmp_path / name).read_bytes()[:4] == b"\x89PNG" for name in names)
