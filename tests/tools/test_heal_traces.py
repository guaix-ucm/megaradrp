import json

import astropy.io.fits as fits
import matplotlib.pyplot
import numpy
import pytest
import yaml

import megaradrp.tools.heal_traces as heal

# heal_traces does not close the file of the traces (#362)
pytestmark = pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning")

INSCONF = "ca3558e3-e50d-4bbc-86bd-da50a0998a48"

# the traces are horizontal, in row 10 + 3 * fibid, 0-based
TRACE_ROWS = {fibid: 10.0 + 3 * fibid for fibid in range(1, 624)}


def trace_state(global_offset=(0.0,)):
    contents = [
        dict(fibid=fibid, boxid=1, start=5, stop=395, fitparms=[TRACE_ROWS[fibid], 0.0]) for fibid in range(1, 624)
    ]
    # fiber 30 without trace
    contents[29]["fitparms"] = []
    return dict(
        instrument="MEGARA",
        tags={"insmode": "LCB", "vph": "LR-R"},
        uuid="11111111-2222-3333-4444-555555555555",
        error_fitting=[30],
        missing_fibers=[],
        total_fibers=623,
        meta_info={"origin": {"insconf_uuid": INSCONF, "date_obs": None}},
        contents=contents,
        type_fqn="megaradrp.products.tracemap.TraceMap",
        boxes_positions=[],
        ref_column=200,
        expected_range=[2, 4092],
        global_offset=list(global_offset),
    )


@pytest.fixture
def files(tmp_path):
    """Image, traces and paths of the outputs of heal_traces"""
    data = numpy.zeros((200, 400), dtype="float32")
    image = tmp_path / "image.fits"
    fits.writeto(image, data)
    traces = tmp_path / "traces.json"
    traces.write_text(json.dumps(trace_state()))
    return tmp_path, image, traces


def run_heal(files, operations=None, extra=()):
    """Run heal_traces with the healing operations, return the updated traces"""
    path, image, traces = files
    args = [str(image), str(traces), "--pdffile", str(path / "out.pdf")]
    args += ["--updated_traces", str(path / "updated.json")]
    if operations is not None:
        healing = path / "healing.yaml"
        healing.write_text(yaml.safe_dump_all(operations))
        args += ["--healing", str(healing)]
    before = matplotlib.pyplot.get_fignums()
    heal.main(args + list(extra))
    assert (path / "out.pdf").exists()
    # the figure is saved without pyplot
    assert matplotlib.pyplot.get_fignums() == before
    return json.loads((path / "updated.json").read_text())


def fitparms(state, fibid):
    return state["contents"][fibid - 1]["fitparms"]


def test_assign_boxes_to_fibers(capsys):
    boxes = [{"name": "1", "nfibers": 2}, {"name": "2b", "nfibers": 3}]
    labels = heal.assign_boxes_to_fibers(boxes, "LCB")
    assert labels == ["1  [1]", "2  [1]", "3  [2b]", "4  [2b]", "5  [2b]"]
    out = capsys.readouterr().out
    assert "Box 2b,  fibers   3 -   5" in out


def test_heal_no_operations(files):
    """Without operations the traces are the same, with a new uuid"""
    state = run_heal(files)
    expected = trace_state()
    assert state["uuid"] != expected["uuid"]
    assert state["contents"] == expected["contents"]


def test_heal_global_offset_argument(files):
    """The global offset of the command line is added to the traces"""
    state = run_heal(files, extra=["--global_offset", "1.5"])
    assert fitparms(state, 1) == [TRACE_ROWS[1] + 1.5, 0.0]
    assert fitparms(state, 30) == []


def test_heal_global_offset_both(files):
    path, image, traces = files
    traces.write_text(json.dumps(trace_state(global_offset=[1.0])))
    with pytest.raises(ValueError, match="global_offset"):
        run_heal(files, extra=["--global_offset", "2.0"])


@pytest.mark.xfail(reason="the global offset of the file is applied twice", strict=True)
def test_heal_global_offset_file(files):
    """The global offset of the file is applied only once"""
    path, image, traces = files
    traces.write_text(json.dumps(trace_state(global_offset=[1.0])))
    state = run_heal(files)
    # the offset is added to the traces, so it must be removed from the file,
    # the extraction applies it again
    assert fitparms(state, 1) == [TRACE_ROWS[1] + 1.0, 0.0]
    assert state["global_offset"] == [0.0]


def test_heal_total_fibers_mismatch(files):
    path, image, traces = files
    state = trace_state()
    state["total_fibers"] = 644
    traces.write_text(json.dumps(state))
    with pytest.raises(ValueError, match="Mismatch"):
        run_heal(files)


def test_heal_does_not_overwrite_traces(files):
    path, image, traces = files
    heal.main([str(image), str(traces), "--pdffile", str(path / "out.pdf"), "--updated_traces", str(traces)])
    assert json.loads(traces.read_text()) == trace_state()


def test_heal_show(files, monkeypatch, capsys):
    """Without PDF file, the plot is shown"""
    path, image, traces = files
    calls = []
    monkeypatch.setattr(heal, "pause_debugplot", lambda *args, **kwargs: calls.append(kwargs))
    heal.main([str(image), str(traces), "--fibids", "--echo"])
    assert calls == [{"pltshow": True, "tight_layout": True}]
    assert "% " in capsys.readouterr().out


@pytest.mark.parametrize(
    "operation, shifted",
    [
        ({"fibid_ini": 3, "fibid_end": 5}, [3, 4, 5]),
        ({"fibid_list": [7, 2]}, [2, 7]),
    ],
)
def test_heal_vertical_shift(files, operation, shifted):
    operation.update(description="vertical_shift_in_pixels", vshift=2.5)
    state = run_heal(files, [operation])
    for fibid in range(1, 10):
        expected = TRACE_ROWS[fibid] + 2.5 if fibid in shifted else TRACE_ROWS[fibid]
        assert fitparms(state, fibid) == [expected, 0.0]


def test_heal_vertical_shift_missing_fiber(files, capsys):
    operation = dict(description="vertical_shift_in_pixels", fibid_list=[30], vshift=2.5)
    state = run_heal(files, [operation], extra=["--verbose"])
    assert fitparms(state, 30) == []
    assert "SKIPPED" in capsys.readouterr().out


@pytest.mark.parametrize(
    "operation, message",
    [
        ({"fibid_list": [1], "vshift": 1, "other": 1}, "Invalid tags"),
        ({"fibid_ini": 5, "fibid_end": 3, "vshift": 1}, "fibid_ini=5 > fibid_end=3"),
        ({"fibid_list": [624], "vshift": 1}, "outside valid range"),
    ],
)
def test_heal_vertical_shift_errors(files, operation, message):
    operation["description"] = "vertical_shift_in_pixels"
    with pytest.raises(ValueError, match=message):
        run_heal(files, [operation])


def test_heal_duplicate_trace(files):
    operation = dict(description="duplicate_trace", fibid_original=4, fibid_duplicated=30, vshift=-1.0)
    state = run_heal(files, [operation], extra=["--verbose"])
    assert fitparms(state, 30) == [TRACE_ROWS[4] - 1.0, 0.0]
    assert fitparms(state, 4) == [TRACE_ROWS[4], 0.0]


def test_heal_duplicate_missing_trace(files, capsys):
    operation = dict(description="duplicate_trace", fibid_original=30, fibid_duplicated=31, vshift=-1.0)
    state = run_heal(files, [operation])
    assert fitparms(state, 31) == [TRACE_ROWS[31], 0.0]
    assert "SKIPPED" in capsys.readouterr().out


@pytest.mark.parametrize("original, duplicated", [(0, 1), (1, 624)])
def test_heal_duplicate_trace_range(files, original, duplicated):
    operation = dict(description="duplicate_trace", fibid_original=original, fibid_duplicated=duplicated, vshift=1)
    with pytest.raises(ValueError, match="outside valid range"):
        run_heal(files, [operation])


@pytest.mark.parametrize("xstart, xstop", [(1, 400), (10, 300)])
def test_heal_extrapolation(files, xstart, xstop):
    operation = dict(description="extrapolation", fibid_ini=1, fibid_end=3, xstart=xstart, xstop=xstop)
    state = run_heal(files, [operation], extra=["--verbose"])
    for fibid in [1, 2, 3]:
        assert state["contents"][fibid - 1]["start"] == xstart
        assert state["contents"][fibid - 1]["stop"] == xstop
    assert state["contents"][3]["start"] == 5


def test_heal_extrapolation_missing_fiber(files, capsys):
    operation = dict(description="extrapolation", fibid_list=[30], xstart=1, xstop=400)
    state = run_heal(files, [operation])
    assert state["contents"][29]["start"] == 5
    assert "SKIPPED" in capsys.readouterr().out


@pytest.mark.parametrize(
    "operation, message",
    [
        ({"fibid_list": [1], "xstart": 1, "xstop": 400, "other": 1}, "Invalid tags"),
        ({"fibid_ini": 5, "fibid_end": 3, "xstart": 1, "xstop": 400}, "fibid_ini=5 > fibid_end=3"),
        ({"fibid_list": [1], "xstart": 400, "xstop": 1}, "xstart=400 > xstop=1"),
        ({"fibid_list": [0], "xstart": 1, "xstop": 400}, "outside valid range"),
    ],
)
def test_heal_extrapolation_errors(files, operation, message):
    operation["description"] = "extrapolation"
    with pytest.raises(ValueError, match=message):
        run_heal(files, [operation])


def test_heal_fit_through_user_points(files):
    """The user points are in image coordinates, starting in 1"""
    operation = dict(
        description="fit_through_user_points",
        fibid=30,
        poldeg=1,
        xstart=10,
        xstop=390,
        user_points=[[11, 101], [201, 111], [391, 120]],
        refit=False,
    )
    state = run_heal(files, [operation], extra=["--verbose"])
    coeff = fitparms(state, 30)
    expected = numpy.polynomial.Polynomial.fit([10, 200, 390], [100, 110, 119], 1).convert().coef
    assert numpy.allclose(coeff, expected)
    assert (state["contents"][29]["start"], state["contents"][29]["stop"]) == (10, 390)


def test_heal_fit_through_user_points_few(files):
    operation = dict(
        description="fit_through_user_points",
        fibid=30,
        poldeg=2,
        xstart=10,
        xstop=390,
        user_points=[[11, 101], [391, 120]],
        refit=False,
    )
    with pytest.raises(ValueError, match="Insufficient number of points"):
        run_heal(files, [operation])


def image_with_ridge(path, row):
    """An image with a gaussian ridge centered in a (0-based) row"""
    yy = numpy.arange(200)[:, numpy.newaxis]
    xx = numpy.arange(400)
    data = 1000 * numpy.exp(-0.5 * ((yy - row(xx)) / 1.2) ** 2)
    fits.writeto(path / "image.fits", data.astype("float32"), overwrite=True)


def test_heal_fit_through_user_points_refit(files, capsys):
    """The trace is refitted to the peaks of the image"""
    path, image, traces = files
    image_with_ridge(path, lambda x: 100.0 + 0 * x)
    operation = dict(
        description="fit_through_user_points",
        fibid=30,
        poldeg=1,
        xstart=10,
        xstop=390,
        user_points=[[11, 103], [391, 103]],
        refit=True,
    )
    state = run_heal(files, [operation])
    # the precision of the refit is checked in test_refit_trace
    coeff = fitparms(state, 30)
    assert abs(coeff[0] - 100.0) < 0.5
    assert abs(coeff[1]) < 1e-6
    assert "refitting...OK!" in capsys.readouterr().out


@pytest.mark.xfail(reason="refit_trace fits the peak in a window centered in the previous row", strict=True)
def test_refit_trace(tmp_path):
    """The peak of a ridge in an integer row is found"""
    data = numpy.zeros((200, 400))
    data[99:102] = [[500], [1000], [500]]
    with fits.HDUList([fits.PrimaryHDU(data)]) as hdul:
        hdul.writeto(tmp_path / "refit.fits")
    coeff = heal.refit_trace(tmp_path / "refit.fits", [100.0, 0.0], 10, 390, 1)
    assert numpy.allclose(coeff, [100.0, 0.0], atol=0.01)


@pytest.mark.xfail(
    reason="refit_trace fits the peaks of the 0-based column j in x = j + 1 (and the window)", strict=True
)
def test_refit_trace_inclined(tmp_path):
    """The refitted trace is in the same coordinates than the traces, x is 0-based"""
    yy = numpy.arange(200)[:, numpy.newaxis]
    xx = numpy.arange(400)
    # a ridge with peak in row 20 + 0.4 * x
    data = numpy.zeros((200, 400))
    for dy, value in [(-1, 500), (0, 1000), (1, 500)]:
        data[yy == numpy.round(20 + 0.4 * xx) + dy] = value
    with fits.HDUList([fits.PrimaryHDU(data)]) as hdul:
        hdul.writeto(tmp_path / "refit.fits")
    coeff = heal.refit_trace(tmp_path / "refit.fits", [20.0, 0.4], 10, 390, 1)
    assert numpy.allclose(coeff, [20.0, 0.4], atol=0.05)


def test_heal_extrapolation_through_user_points(files):
    operation = dict(
        description="extrapolation_through_user_points",
        fibid=2,
        xstart_reuse=100,
        xstop_reuse=300,
        resampling=20,
        poldeg=1,
        xstart=1,
        xstop=400,
        user_points=[[396, TRACE_ROWS[2] + 1]],
        refit=False,
    )
    state = run_heal(files, [operation], extra=["--verbose", "--fibids"])
    coeff = fitparms(state, 2)
    assert numpy.allclose(coeff, [TRACE_ROWS[2], 0.0], atol=1e-6)
    assert (state["contents"][1]["start"], state["contents"][1]["stop"]) == (1, 400)


def test_heal_extrapolation_through_user_points_refit(files):
    path, image, traces = files
    image_with_ridge(path, lambda x: TRACE_ROWS[2] - 2 + 0 * x)
    operation = dict(
        description="extrapolation_through_user_points",
        fibid=2,
        xstart_reuse=100,
        xstop_reuse=300,
        resampling=20,
        poldeg=1,
        xstart=10,
        xstop=390,
        user_points=[[391, TRACE_ROWS[2] + 1]],
        refit=True,
    )
    state = run_heal(files, [operation])
    coeff = fitparms(state, 2)
    assert abs(coeff[0] - (TRACE_ROWS[2] - 2)) < 0.5
    assert abs(coeff[1]) < 1e-6


@pytest.mark.parametrize(
    "operation, message",
    [
        ({"xstart_reuse": 300, "xstop_reuse": 100, "xstart": 1, "xstop": 400}, "xstart_reuse=300 > xstop_reuse=100"),
        ({"xstart_reuse": 100, "xstop_reuse": 300, "xstart": 400, "xstop": 1}, "xstart=400 > xstop=1"),
    ],
)
def test_heal_extrapolation_through_user_points_errors(files, operation, message):
    operation.update(
        description="extrapolation_through_user_points", fibid=2, resampling=20, poldeg=1, user_points=[], refit=False
    )
    with pytest.raises(ValueError, match=message):
        run_heal(files, [operation])


def test_heal_sandwich(files):
    """The trace is interpolated between the neighbours, the fiber is removed from error_fitting"""
    operation = dict(description="sandwich", fibid=30, fraction=0.25, neighbours=[29, 31], xstart=5, xstop=395)
    state = run_heal(files, [operation], extra=["--verbose"])
    assert fitparms(state, 30) == [TRACE_ROWS[29] + 0.25 * (TRACE_ROWS[31] - TRACE_ROWS[29]), 0.0]
    assert state["error_fitting"] == []


@pytest.mark.parametrize(
    "operation, message",
    [
        ({"neighbours": [29, 31], "xstart": 395, "xstop": 5}, "xstart=395 > xstop=5"),
    ],
)
def test_heal_sandwich_errors(files, operation, message):
    operation.update(description="sandwich", fibid=30, fraction=0.5)
    with pytest.raises(ValueError, match=message):
        run_heal(files, [operation])


def test_heal_sandwich_neighbours(files):
    path, image, traces = files
    state = trace_state()
    state["contents"][28]["fibid"] = 128
    traces.write_text(json.dumps(state))
    operation = dict(description="sandwich", fibid=30, fraction=0.5, neighbours=[29, 31], xstart=5, xstop=395)
    with pytest.raises(ValueError, match="Unexpected fiber numbers"):
        run_heal(files, [operation])


def test_heal_renumber_down(files):
    """With shift -1 the fibers move to the previous fibid, the last one is empty"""
    operation = dict(description="renumber_fibids_within_box", fibid_ini=23, fibid_end=25, fibid_shift=-1)
    state = run_heal(files, [operation], extra=["--verbose"])
    for fibid in [22, 23, 24]:
        assert state["contents"][fibid - 1]["fibid"] == fibid
        assert fitparms(state, fibid) == [TRACE_ROWS[fibid + 1], 0.0]
    assert fitparms(state, 25) == []


def test_heal_renumber_up(files):
    """With shift +1 the fibers move to the next fibid, the first one is empty"""
    operation = dict(description="renumber_fibids_within_box", fibid_ini=23, fibid_end=25, fibid_shift=1)
    state = run_heal(files, [operation])
    for fibid in [24, 25, 26]:
        assert state["contents"][fibid - 1]["fibid"] == fibid
        assert fitparms(state, fibid) == [TRACE_ROWS[fibid - 1], 0.0]
    assert fitparms(state, 23) == []


@pytest.mark.parametrize(
    "operation, message",
    [
        ({"fibid_ini": 25, "fibid_end": 23, "fibid_shift": 1}, "fibid_ini=25 > fibid_end=23"),
        # fibers 1-21 are in box 8b, 22-42 in box 7b
        ({"fibid_ini": 20, "fibid_end": 23, "fibid_shift": 1}, "different fiber boxes"),
        ({"fibid_ini": 23, "fibid_end": 25, "fibid_shift": 2}, "must be -1 or 1"),
    ],
)
def test_heal_renumber_errors(files, operation, message):
    operation["description"] = "renumber_fibids_within_box"
    with pytest.raises(ValueError, match=message):
        run_heal(files, [operation])


@pytest.mark.xfail(reason="the box is compared with label[4:], that depends on the digits of fibid", strict=True)
def test_heal_renumber_digits(files):
    """Fibers 5 to 15 are in the same box"""
    operation = dict(description="renumber_fibids_within_box", fibid_ini=5, fibid_end=15, fibid_shift=1)
    state = run_heal(files, [operation])
    assert fitparms(state, 5) == []


def test_heal_several_operations(files):
    """The operations are applied in order"""
    operations = [
        dict(description="vertical_shift_in_pixels", fibid_list=[2], vshift=1.0),
        dict(description="duplicate_trace", fibid_original=2, fibid_duplicated=30, vshift=1.0),
    ]
    state = run_heal(files, operations)
    assert fitparms(state, 30) == [TRACE_ROWS[2] + 2.0, 0.0]


@pytest.mark.parametrize(
    "operation, error, message",
    [
        ({"fibid": 1}, ValueError, "'description' not found"),
        ({"description": "other"}, ValueError, "Unexpected healing method: other"),
        ({"description": "sandwich", "fibid": 1}, SystemExit, "Missing tags found in 'sandwich'"),
        (
            {"description": "duplicate_trace", "fibid_original": 1, "fibid_duplicated": 2, "vshift": 1, "other": 1},
            SystemExit,
            "Unexpected tags found",
        ),
        (
            {"description": "duplicate_trace", "fibid_original": 1, "vshift": 1, "other": 1},
            SystemExit,
            "Unexpected and missing tags found",
        ),
    ],
)
def test_heal_wrong_operations(files, operation, error, message):
    with pytest.raises(error, match=message):
        run_heal(files, [operation])
