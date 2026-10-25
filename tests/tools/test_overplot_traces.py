import json

import astropy.io.fits as fits
import numpy
import pytest
from numina.util.jsonencoder import ExtEncoder

import megaradrp.tools.overplot_traces as overplot
from megaradrp.testing.create_tracemap import create_test_tracemap2

INSCONF = "ca3558e3-e50d-4bbc-86bd-da50a0998a48"


def test_overplot_traces_help(capsys):
    """Check that the program runs with --help"""
    with pytest.raises(SystemExit):
        overplot.main(["--help"])

    out, err = capsys.readouterr()
    assert "overplot traces" in out


@pytest.fixture
def files(tmp_path):
    """Image and traces, fiber 30 without trace"""
    fits.writeto(tmp_path / "image.fits", numpy.zeros((200, 400), dtype="float32"))
    tracemap = create_test_tracemap2()
    tracemap.meta_info["origin"] = {"insconf_uuid": INSCONF, "date_obs": None}
    tracemap.contents[29].fitparms = []
    tracemap.global_offset = numpy.polynomial.Polynomial([0.5])
    with open(tmp_path / "traces.json", "w") as fd:
        json.dump(tracemap.__getstate__(), fd, cls=ExtEncoder)
    return tmp_path


def run_overplot(path, extra=()):
    args = [str(path / "image.fits"), str(path / "traces.json"), "--pdffile", str(path / "out.pdf")]
    overplot.main(args + list(extra))
    assert (path / "out.pdf").exists()


def plotted_traces(monkeypatch):
    """Record the arguments of the traces plotted"""
    calls = []
    plot_trace = overplot.plot_trace

    def recorder(ax, center_model, *args, **kwargs):
        calls.append((center_model, args, kwargs))
        plot_trace(ax, center_model, *args, **kwargs)

    monkeypatch.setattr(overplot, "plot_trace", recorder)
    return calls


def test_overplot_traces(files, monkeypatch, capsys):
    calls = plotted_traces(monkeypatch)
    run_overplot(files, ["--fibids"])
    # fiber 30 is missing
    assert len(calls) == 622
    assert "Missing fiber: 30  [7b]" in capsys.readouterr().out
    center_model, args, kwargs = calls[0]
    # start, stop, ix_offset, rawimage, fibids, fiblabel
    assert args == (4, 4090, 0, False, True, "1  [8b]")
    # the global offset of the file
    assert kwargs["correction"] == 0.5
    assert center_model(100) == 203.5


def test_overplot_traces_options(files, monkeypatch):
    calls = plotted_traces(monkeypatch)
    run_overplot(files, ["--rawimage", "--global-offset", "1", "0.01"])
    center_model, args, kwargs = calls[0]
    assert args[2:4] == (50, True)
    # the offset depends on the position of the trace in the reference column
    assert kwargs["correction"] == pytest.approx(1 + 0.01 * 203.5)


def test_overplot_traces_show(files, monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(overplot, "pause_debugplot", lambda *args, **kwargs: calls.append(kwargs))
    overplot.main([str(files / "image.fits"), str(files / "traces.json"), "--echo"])
    assert calls == [{"pltshow": True, "tight_layout": True}]
    assert "% " in capsys.readouterr().out


def test_overplot_traces_malformed(files):
    state = json.loads((files / "traces.json").read_text())
    del state["type_fqn"]
    (files / "traces.json").write_text(json.dumps(state))
    with pytest.raises(ValueError, match="'type_fqn' missing"):
        run_overplot(files)


def test_overplot_traces_total_fibers(files):
    state = json.loads((files / "traces.json").read_text())
    state["total_fibers"] = 644
    (files / "traces.json").write_text(json.dumps(state))
    with pytest.raises(ValueError, match="Mismatch"):
        run_overplot(files)


@pytest.mark.parametrize("rawimage, expected", [(False, [2000.5, 2100.5]), (True, [2000.5, 2200.5])])
def test_plot_trace_rawimage(rawimage, expected):
    """In raw images, the traces in the upper half are moved 100 pixels up"""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    model = numpy.polynomial.Polynomial([1999.0, 0.5])
    overplot.plot_trace(ax, model, 1, 201, 0, rawimage, False, "", "blue")
    (line,) = ax.get_lines()
    plt.close(fig)
    xdata, ydata = line.get_data()
    assert (xdata[0], xdata[-1]) == (1, 201)
    assert [ydata[0], ydata[-1]] == expected


def test_plot_trace_full_range():
    """With start and stop 0, the trace is drawn in the whole detector"""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.set_ylim(0, 300)
    overplot.plot_trace(ax, numpy.polynomial.Polynomial([100.0]), 0, 0, 0, False, True, "label", "blue")
    (line,) = ax.get_lines()
    (text,) = ax.texts
    plt.close(fig)
    assert len(line.get_xdata()) == 4096
    assert text.get_position() == (2048, 100.0)
