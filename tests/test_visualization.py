import astropy.io.fits as fits
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy
import pytest
from astropy.wcs import WCS

import megaradrp.visualization as vis
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.wavecalibration import rss_add_wcs
from megaradrp.testing.create_image import create_rss


@pytest.fixture(autouse=True)
def close_figures(monkeypatch):
    """The figures are not shown, they are closed after the test"""
    monkeypatch.setattr(plt, "show", lambda: None)
    yield
    plt.close("all")


def test_hexplot():
    fig, ax = plt.subplots()
    x = numpy.array([0.0, 1.0, 2.0])
    y = numpy.array([0.0, 0.5, 0.0])
    z = numpy.array([1.0, 2.0, 3.0])
    col = vis.hexplot(ax, x, y, z, scale=0.5, vmin=0, vmax=4)
    assert isinstance(col, mcoll.PolyCollection)
    assert col in ax.collections
    assert numpy.all(col.get_offsets() == numpy.column_stack([x, y]))
    assert numpy.all(col.get_array() == z)
    assert col.get_clim() == (0, 4)
    # an hexagon with flat top and bottom
    (polygon,) = col.get_paths()
    assert polygon.vertices[:, 1].max() == pytest.approx(0.99 * 0.5 * 0.5)
    # the limits contain the hexagons
    xmin, xmax = ax.get_xlim()
    assert xmin < 0 - 0.5 and xmax > 2 + 0.5


def test_hexplot_lognorm():
    """With logarithmic normalization the zeros are avoided"""
    fig, ax = plt.subplots()
    z = numpy.array([0.0, 9.0, 99.0])
    col = vis.hexplot(ax, [0, 1, 2], [0, 0, 0], z, norm=mcolors.LogNorm())
    assert numpy.all(col.get_array() == [1.0, 10.0, 100.0])
    assert col.norm.vmin == 1.0
    assert col.norm.vmax == 100.0


def test_hexplot_extent():
    fig, ax = plt.subplots()
    vis.hexplot(ax, [0, 1], [0, 0], numpy.array([1.0, 2.0]), extent=(-10, 10, -5, 5))
    xmin, xmax = ax.get_xlim()
    assert xmin <= -10
    assert xmax >= 10


def test_is_inside():
    assert vis.is_inside_1d(10, [0, 9])
    assert not vis.is_inside_1d(10, [0, 10])
    assert not vis.is_inside_1d(10, [-1])
    assert vis.is_inside((10, 20), [[0, 9], [0, 19]])
    assert not vis.is_inside((10, 20), [[0, 9], [0, 20]])


def wcs_wl(coordinate_type, header="AWAV"):
    """The WCS of the RSS of MEGARA (header "AWAV"), without CUNIT1, or before megaradrp 0.9 (CTYPE1 WAVELENGTH)"""
    if coordinate_type == "wcs":
        hdr = rss_add_wcs(fits.Header(), 4000.0, 0.5, 1.0)
        if header == "no-cunit":
            del hdr["CUNIT1"]
        elif header == "WAVELENGTH":
            hdr["CTYPE1"] = "WAVELENGTH"
        return vis.spectral_wcs(hdr)
    return WCS(naxis=2)


@pytest.fixture
def rssdata():
    """Each fiber has its row number, the column number is added"""
    return numpy.arange(3)[:, numpy.newaxis] * 1000.0 + numpy.arange(100)


@pytest.mark.parametrize(
    "coordinate_type, column, region",
    [("pixel", 21, [11, 31]), ("wcs", 4010.0, [4005.0, 4015.0])],
)
def test_extract_zval(rssdata, coordinate_type, column, region):
    """The pixel coordinates start in 1"""
    wcs = wcs_wl(coordinate_type)
    zval = vis.extract_zval(rssdata, wcs, coordinate_type, column=column, average_region=region)
    assert numpy.all(zval == [20, 1020, 2020])
    zval = vis.extract_zval(rssdata, wcs, coordinate_type, average_region=region)
    # columns 10 to 30, the limits are included
    assert numpy.all(zval == [20, 1020, 2020])


@pytest.mark.parametrize("header", ["AWAV", "no-cunit", "WAVELENGTH"])
def test_extract_zval_angstrom(rssdata, header):
    """The wavelengths are in Angstrom, astropy converts those of the WCS to meters (#363)"""
    wcs = wcs_wl("wcs", header)
    zval = vis.extract_zval(rssdata, wcs, "wcs", column=4010.0)
    assert numpy.all(zval == [20, 1020, 2020])
    zval = vis.extract_zval(rssdata, wcs, "wcs", average_region=[4005.0, 4015.0])
    assert numpy.all(zval == [20, 1020, 2020])


def test_extract_zval_unit(rssdata):
    wcs = wcs_wl("wcs")
    zval = vis.extract_zval(rssdata, wcs, "wcs", column=401.0, unit="nm")
    assert numpy.all(zval == [20, 1020, 2020])


def test_extract_zval_column(rssdata):
    """The column is enough, without region (#363)"""
    zval = vis.extract_zval(rssdata, wcs_wl("pixel"), "pixel", column=21)
    assert numpy.all(zval == [20, 1020, 2020])


@pytest.mark.parametrize("wl, pixel", [(4010.0, 20), (4010.2, 20), (4010.3, 21), (4009.8, 20)])
def test_transform_to_pixel_nearest(wl, pixel):
    """The index of the nearest pixel"""
    assert list(vis.transform_to_pixel(wcs_wl("wcs"), [wl], unit="Angstrom")) == [pixel]


def test_extract_zval_continuum(rssdata):
    wcs = wcs_wl("pixel")
    zval = vis.extract_zval(rssdata, wcs, "pixel", average_region=[11, 31], continuum_region=[51, 71])
    assert numpy.all(zval == -40)


def test_extract_zval_1d():
    data = numpy.array([1.0, 2.0])
    assert vis.extract_zval(data, None, "pixel", average_region=[1, 2]) is data


@pytest.mark.parametrize(
    "data, kwargs, message",
    [
        (numpy.zeros((3, 100)), {"continuum_region": [1, 2]}, "are None"),
        (numpy.zeros((2, 3, 100)), {"average_region": [1, 2]}, "1D or 2D"),
        (numpy.zeros((3, 100)), {"average_region": [11, 200]}, "region to average is outside image"),
        (numpy.zeros((3, 100)), {"column": 101, "average_region": [1, 2]}, "column is outside image"),
    ],
)
def test_extract_zval_errors(data, kwargs, message):
    with pytest.raises(ValueError, match=message):
        vis.extract_zval(data, wcs_wl("pixel"), "pixel", **kwargs)


def lcb_rss(path, name="rss.fits"):
    """A LCB RSS, the value of each fiber is its fibid"""
    data = numpy.repeat(numpy.arange(1, 624, dtype="float32")[:, numpy.newaxis], 4000, axis=1)
    rss = create_rss(data, numpy.ones_like(data))
    rss_add_wcs(rss[0].header, 4000.0, 0.5, 1.0)
    rss[0].header["IPA"] = 90.0
    fname = path / name
    rss.writeto(fname)
    return fname


def displayed_values():
    """The values of the hexagons of the current figure"""
    ax = plt.gcf().axes[0]
    (col,) = [col for col in ax.collections if isinstance(col, mcoll.PolyCollection)]
    return col.get_array()


@pytest.fixture
def fp_conf():
    return FocalPlaneConf.from_img(create_rss(numpy.zeros((623, 1)), numpy.zeros((623, 1))))


def test_main(tmp_path, fp_conf):
    """The sky fibers are not displayed"""
    fname = lcb_rss(tmp_path)
    vis.main([str(fname), "--title", "title", "--label", "label"])
    skyfibers = set(fp_conf.sky_fibers())
    expected = [fibid for fibid in range(1, 624) if fibid not in skyfibers]
    assert numpy.all(displayed_values() == expected)
    ax = plt.gcf().axes[0]
    assert ax.get_title() == "title"


def test_main_options(tmp_path):
    fname = lcb_rss(tmp_path)
    args = ["--plot-sky", "--colormap", "viridis", "--stretch", "log", "--hex-rel-size", "0.5"]
    args += ["--column", "100", "--plot-nominal-config"]
    vis.main(args + [str(fname)])
    assert numpy.all(displayed_values() == numpy.arange(1, 624))


def test_main_hide_values(tmp_path):
    fname = lcb_rss(tmp_path)
    vis.main(["--plot-sky", "--hide-values", "--min-cut", "100", "--max-cut", "200", str(fname)])
    # the hexagons out of the range are not drawn
    assert numpy.all(displayed_values() == numpy.arange(100, 201))


def test_main_fix_missing(tmp_path):
    """Fiber 623 is interpolated"""
    fname = lcb_rss(tmp_path)
    vis.main(["--plot-sky", "--fix-missing", str(fname)])
    assert displayed_values()[-1] != 623


# the FIBERS extension has no data
@pytest.mark.filterwarnings("ignore:The WCS transformation has more axes:astropy.wcs.wcs.FITSFixedWarning")
@pytest.mark.parametrize("pa_from_header", [False, True])
def test_main_wcs_grid(tmp_path, pa_from_header):
    fname = lcb_rss(tmp_path)
    args = ["--wcs-grid", str(fname)]
    if pa_from_header:
        args.append("--wcs-pa-from-header")
    vis.main(args)
    ax = plt.gcf().axes[0]
    assert hasattr(ax, "coords")


def test_main_contour(tmp_path, capsys):
    fname = lcb_rss(tmp_path)
    save = tmp_path / "contour.fits"
    vis.main(["--contour", "--contour-levels", "[100, 300]", "--contour-image-save", str(save), str(fname)])
    ax = plt.gcf().axes[0]
    assert "contour levels [100. 300.]" in capsys.readouterr().out
    with fits.open(save) as hdul:
        assert hdul[0].data.shape == (623, 1)
    assert len(ax.collections) > 1


def test_main_contour_image(tmp_path):
    fname = lcb_rss(tmp_path)
    contour = lcb_rss(tmp_path, "contour.fits")
    args = ["--contour", "--contour-image", str(contour), "--contour-image-column", "100", "--contour-is-density"]
    vis.main(args + [str(fname)])
    ax = plt.gcf().axes[0]
    assert len(ax.collections) > 1


def test_main_display_fibid(tmp_path, monkeypatch):
    calls = []
    text_path = vis.TextPath

    def recorder(xy, s, size):
        calls.append(s)
        return text_path(xy, s, size=size)

    monkeypatch.setattr(vis, "TextPath", recorder)
    fname = lcb_rss(tmp_path)
    vis.main(["--plot-sky", "--display-fibid", str(fname)])
    # number and name of each spaxel
    assert len(calls) == 2 * 623


def test_main_display_fibid_labels(tmp_path, monkeypatch, fp_conf):
    """Each spaxel has the label of its fibid and its name (#363)"""
    calls = []
    text_path = vis.TextPath

    def recorder(xy, s, size):
        calls.append(s)
        return text_path(xy, s, size=size)

    monkeypatch.setattr(vis, "TextPath", recorder)
    fname = lcb_rss(tmp_path)
    vis.main(["--display-fibid", str(fname)])
    skyfibers = set(fp_conf.sky_fibers())
    expected = [f"{fibid:03d}" for fibid in range(1, 624) if fibid not in skyfibers]
    assert calls[::2] == expected
    names = {fiber.fibid: fiber.name for fiber in fp_conf.fibers.values()}
    assert calls[1::2] == [names[int(fibid)] for fibid in expected]


def test_main_wcs_coordinates(tmp_path):
    """The regions in Angstrom, the RSS has CUNIT1 (#363)"""
    fname = lcb_rss(tmp_path)
    vis.main(["--plot-sky", "--coordinate-type", "wcs", "--average-region", "4500", "4510.5", str(fname)])
    assert numpy.all(displayed_values() == numpy.arange(1, 624))
    vis.main(["--plot-sky", "--coordinate-type", "wcs", "--column", "5999.5", str(fname)])
    assert numpy.all(displayed_values() == numpy.arange(1, 624))


@pytest.mark.parametrize("args", [[], ["--contour", "--contour-image", "contour.fits", "--column", "4500"]])
def test_main_wcs_coordinates_required(tmp_path, args, capsys):
    """The default regions are in pixels, with wcs the region is required"""
    fname = lcb_rss(tmp_path)
    with pytest.raises(SystemExit):
        vis.main(["--coordinate-type", "wcs"] + args + [str(fname)])
    assert "--coordinate-type wcs requires" in capsys.readouterr().err
