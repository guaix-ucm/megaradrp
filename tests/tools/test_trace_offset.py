import json

import astropy.io.fits as fits
import numpy
import pytest
from numina.util.jsonencoder import ExtEncoder

import megaradrp.tools.trace_offset as trace_offset
from megaradrp.testing.create_traceoffset import create_offset_image, create_offset_tracemap

OPTS = ["--nbands", "4", "--band-width", "32"]


def write_files(path, flux=None):
    tracemap = create_offset_tracemap()
    image = create_offset_image(tracemap, -2.4, flux=flux)
    header = fits.Header({"BUNIT": "ELECTRON", "RDNOISE1": 3.4})
    fits.writeto(path / "image.fits", image.astype("float32"), header)
    with open(path / "traces.json", "w") as fd:
        json.dump(tracemap.__getstate__(), fd, cls=ExtEncoder)
    return [str(path / "image.fits"), str(path / "traces.json")]


def test_trace_offset_help(capsys):
    with pytest.raises(SystemExit):
        trace_offset.main(["--help"])
    assert "extraction_offset" in capsys.readouterr().out


def test_trace_offset(tmp_path, capsys):
    args = write_files(tmp_path)
    plot, out = tmp_path / "offset.png", tmp_path / "offset.json"
    status = trace_offset.main(args + OPTS + ["--plot", str(plot), "--json", str(out)])
    assert status == 0
    text = capsys.readouterr().out
    assert "extraction_offset: [-2.4" in text
    assert "reliable: True" in text
    assert plot.exists()
    result = json.loads(out.read_text())
    assert result["coef"][0] == pytest.approx(-2.4, abs=0.03)
    assert result["reliable"]
    assert len(result["bands"]) == 4


def test_trace_offset_not_reliable(tmp_path, capsys):
    """Without signal, the result is not reliable and the exit status is 1"""
    args = write_files(tmp_path, flux=numpy.zeros(40))
    assert trace_offset.main(args + OPTS + ["--no-slope", "--ron", "3.4"]) == 1
    assert "reliable: False" in capsys.readouterr().out
