import copy
import math
import warnings

import astropy.units as u
import astropy.wcs
from astropy.io import fits
import numpy as np
import pytest
from scipy.spatial import cKDTree

import megaradrp.instrument.constants as cons
from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.processing.cube import (
    HEX_SCALE,
    calc_matrix_from_fiberconf,
    create_cube,
    create_cube_from_array,
    create_cube_from_rss,
    main,
    merge_wcs,
)
from megaradrp.processing.hexgrid import HA_HEX, hexgrid_extremes
from megaradrp.processing.wavecalibration import header_add_barycentric_correction
from megaradrp.processing.wcs import compute_pa_from_ipa
from megaradrp.testing.create_header import create_sky_header2, create_spec_header2


@pytest.fixture(scope="module")
def lcb_fiber_header():
    """Header of FIBERS of the LCB, with a sky WCS without rotation"""
    hdrf = create_default_fiber_header("LCB")
    sky = create_sky_header2()
    for key in sky:
        if key != "DATE-OBS":
            hdrf[key] = sky[key]
    hdrf["PC1_1"], hdrf["PC1_2"], hdrf["PC2_1"], hdrf["PC2_2"] = 1.0, 0.0, 0.0, 1.0
    return hdrf


@pytest.fixture
def fiber_header(lcb_fiber_header):
    return lcb_fiber_header.copy()


@pytest.fixture(scope="module")
def focal_plane(lcb_fiber_header):
    return FocalPlaneConf.from_header(lcb_fiber_header)


def make_rss(data, fiber_header):
    """A RSS of the LCB with data, the spectral WCS of create_spec_header2"""
    prim = fits.PrimaryHDU(np.asarray(data, dtype="float32"), header=create_spec_header2())
    return fits.HDUList([prim, fits.ImageHDU(header=fiber_header, name="FIBERS")])


def connected_rows(focal_plane):
    return [fiber.fibid - 1 for fiber in focal_plane.connected_fibers()]


def test_calc_matrix_from_fiberconf(focal_plane):
    r0l, _ = calc_matrix_from_fiberconf(focal_plane)
    fibids = [fiber.fibid for fiber in focal_plane.connected_fibers()]
    assert r0l.shape == (2, len(fibids))
    # the reference fiber is the origin
    assert np.allclose(r0l[:, fibids.index(614)], 0.0)
    # the fibers are in a hexagonal grid of unit spacing
    distances, _ = cKDTree(r0l.T).query(r0l.T, k=2)
    assert np.allclose(distances[:, 1], 1.0, atol=2e-3)


def test_calc_matrix_from_fiberconf_units(focal_plane):
    """The grid is the same with the positions of the fibers in mm or in arcsec"""
    r0l_mm, ref_mm = calc_matrix_from_fiberconf(focal_plane)
    fp_arcsec = copy.deepcopy(focal_plane)
    fp_arcsec.funit = "arcsec"
    factor = HEX_SCALE / cons.SPAXEL_SCALE.to(u.mm).value
    for fiber in fp_arcsec.fibers.values():
        fiber.x *= factor
        fiber.y *= factor
    r0l_arcsec, ref_arcsec = calc_matrix_from_fiberconf(fp_arcsec)
    assert np.allclose(r0l_mm, r0l_arcsec)
    assert np.allclose(ref_mm, ref_arcsec)


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)


@pytest.mark.parametrize("p", [1, 2])
def test_create_cube_from_array_shape(focal_plane, p):
    nwl = 3
    scale = 1.0
    cube = create_cube_from_array(np.ones((623, nwl)), focal_plane, p=p, target_scale_arcsec=scale)
    r0l, _ = calc_matrix_from_fiberconf(focal_plane)
    (i1min, i1max), (j1min, j1max) = hexgrid_extremes(r0l, scale / HEX_SCALE)
    # the wavelength is the first axis, the last in FITS
    assert cube.shape == (nwl, i1max - i1min + 1, j1max - j1min + 1)
    assert cube.dtype == np.float32


def test_create_cube_from_array_1d(focal_plane):
    rss = np.arange(623.0)
    cube = create_cube_from_array(rss, focal_plane)
    cube2d = create_cube_from_array(rss[:, np.newaxis], focal_plane)
    assert cube.shape[0] == 1
    assert np.allclose(cube, cube2d)


@pytest.mark.parametrize("p", [1, 2])
def test_create_cube_from_array_conserve_flux(focal_plane, p):
    """The flux is multiplied by the area of the pixel over the area of the fiber"""
    scale = 1.0
    rss = np.ones((623, 1))
    cube = create_cube_from_array(rss, focal_plane, p=p, target_scale_arcsec=scale)
    cube_nc = create_cube_from_array(rss, focal_plane, p=p, target_scale_arcsec=scale, conserve_flux=False)
    target_scale = scale / HEX_SCALE
    assert np.allclose(cube, cube_nc * target_scale**2 / HA_HEX, rtol=1e-5)


@pytest.mark.parametrize("scale", [0.3, 0.5, 1.0])
@pytest.mark.parametrize("p", [1, 2])
def test_create_cube_from_array_total_flux(focal_plane, p, scale):
    """The flux of the cube is the flux of the fibers (#343)"""
    rss = np.ones((623, 1))
    cube = create_cube_from_array(rss, focal_plane, p=p, target_scale_arcsec=scale)
    total = rss[connected_rows(focal_plane)].sum()
    assert math.isclose(cube.sum(), total, rel_tol=0.02)


@pytest.mark.parametrize("p", [1, 2])
@pytest.mark.parametrize("fibid", [1, 311, 614])
def test_create_cube_from_array_fiber_position(focal_plane, p, fibid):
    """The light of a fiber is in the pixel nearest to its position"""
    scale = 0.5
    rss = np.zeros(623)
    rss[fibid - 1] = 1.0
    img = create_cube_from_array(rss, focal_plane, p=p, target_scale_arcsec=scale)[0]

    r0l, _ = calc_matrix_from_fiberconf(focal_plane)
    fibids = [fiber.fibid for fiber in focal_plane.connected_fibers()]
    x, y = r0l[:, fibids.index(fibid)]
    target_scale = scale / HEX_SCALE
    (i1min, _), (j1min, _) = hexgrid_extremes(r0l, target_scale)
    expected = (round(y / target_scale) - i1min, round(x / target_scale) - j1min)
    assert np.unravel_index(np.argmax(img), img.shape) == expected


def test_create_cube_from_rss_header(fiber_header):
    nwl = 5
    scale = 1.0
    rss = make_rss(np.ones((623, nwl)), fiber_header)
    cube = create_cube_from_rss(rss, target_scale_arcsec=scale)
    hdr = cube[0].header
    hdr_spec = create_spec_header2()

    assert cube[0].data.shape[0] == nwl
    # the extensions are copied
    assert "FIBERS" in cube
    assert [hdr["CTYPE1"], hdr["CTYPE2"], hdr["CTYPE3"]] == ["RA---TAN", "DEC--TAN", hdr_spec["CTYPE1"]]
    # sky axes from FIBERS, with the size of the pixel of the cube
    assert np.allclose([hdr["CRVAL1"], hdr["CRVAL2"]], [fiber_header["CRVAL1"], fiber_header["CRVAL2"]])
    assert np.allclose([hdr["CDELT1"], hdr["CDELT2"]], [-scale / 3600, scale / 3600])
    # spectral axis from the primary header
    assert np.allclose([hdr["CRVAL3"], hdr["CDELT3"]], [hdr_spec["CRVAL1"], hdr_spec["CDELT1"]])


@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("scale", [0.3, 0.7])
def test_create_cube_from_rss_wcs(fiber_header, focal_plane, rotated, scale):
    """The WCS of the cube gives the sky position of the fibers (#342)"""
    if rotated:
        sky = create_sky_header2()
        for key in ["PC1_1", "PC1_2", "PC2_1", "PC2_2"]:
            fiber_header[key] = sky[key]
    rss = make_rss(np.ones((623, 1)), fiber_header)
    cube = create_cube_from_rss(rss, target_scale_arcsec=scale)
    wcs_cube = astropy.wcs.WCS(cube[0].header, fix=False).celestial
    wcs_fibers = astropy.wcs.WCS(fiber_header, fix=False)

    r0l, _ = calc_matrix_from_fiberconf(focal_plane)
    fibids = [fiber.fibid for fiber in focal_plane.connected_fibers()]
    target_scale = scale / HEX_SCALE
    (i1min, _), (j1min, _) = hexgrid_extremes(r0l, target_scale)
    for fibid in [1, 311, 614]:
        fiber = focal_plane.fibers[fibid]
        # the positions of the fibers are the pixel coordinates of FIBERS
        expected = wcs_fibers.all_pix2world(fiber.x, fiber.y, 1)
        # position of the fiber in the array of the cube, 0-based
        x, y = r0l[:, fibids.index(fibid)]
        world = wcs_cube.all_pix2world(x / target_scale - j1min, y / target_scale - i1min, 0)
        assert np.allclose(world, expected, rtol=0, atol=1e-3 / 3600)


def test_sub_wcs():
    hdr_sky = create_sky_header2()
    hdr_spec = create_spec_header2()
    hdr_spec = header_add_barycentric_correction(hdr_spec)
    wcs_sky = astropy.wcs.WCS(header=hdr_sky, fix=False)
    wcs_spec = astropy.wcs.WCS(header=hdr_spec, key="B", fix=False)
    wcs3 = wcs_sky.sub([1, 2, 0])
    wcs3.wcs.ctype[2] = wcs_spec.wcs.ctype[0]
    wcs3.wcs.crval[2] = wcs_spec.wcs.crval[0]
    wcs3.wcs.crpix[2] = wcs_spec.wcs.crpix[0]
    wcs3.wcs.cdelt[2] = wcs_spec.wcs.cdelt[0]
    wcs3.wcs.cunit[2] = wcs_spec.wcs.cunit[0]
    wcs3.wcs.specsys = wcs_spec.wcs.specsys
    wcs3.wcs.ssysobs = wcs_spec.wcs.ssysobs
    wcs3.wcs.velosys = wcs_spec.wcs.velosys
    # hdr3 = wcs3.to_header(key='B')
    assert True


@pytest.mark.xfail(
    reason="Astropy issue #15814: spurious END warning from find_all_wcs()",
    strict=False,
)
def test_find_all_wcs_should_not_warn():
    hdr_spec1 = create_spec_header2()

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        astropy.wcs.find_all_wcs(hdr_spec1)

    assert not w


def test_merge_wcs():
    hdr_spec = create_spec_header2()
    hdr_spec = header_add_barycentric_correction(hdr_spec)
    hdr_sky = create_sky_header2()
    out = hdr_spec.copy()
    merge_wcs(hdr_sky, hdr_spec, out=out)

    wcs0 = astropy.wcs.WCS(header=out, fix=False)
    wcsB = astropy.wcs.WCS(header=out, key="B", fix=False)

    assert list(wcs0.wcs.ctype) == ["RA---TAN", "DEC--TAN", "AWAV"]
    assert np.allclose(
        [out["CRVAL1"], out["CRVAL2"], out["CRVAL3"]],
        [hdr_sky["CRVAL1"], hdr_sky["CRVAL2"], hdr_spec["CRVAL1"]],
    )
    assert np.allclose(
        [out["CDELT1"], out["CDELT2"], out["CDELT3"]],
        [hdr_sky["CDELT1"], hdr_sky["CDELT2"], hdr_spec["CDELT1"]],
    )
    assert list(wcsB.wcs.ctype) == ["RA---TAN", "DEC--TAN", "AWAV"]
    assert np.allclose(
        [out["CRVAL1B"], out["CRVAL2B"], out["CRVAL3B"]],
        [hdr_sky["CRVAL1"], hdr_sky["CRVAL2"], hdr_spec["CRVAL1B"]],
    )
    assert np.allclose(
        [out["CDELT1B"], out["CDELT2B"], out["CDELT3B"]],
        [hdr_sky["CDELT1"], hdr_sky["CDELT2"], hdr_spec["CDELT1B"]],
    )


def test_merge_wcs_only_primary():
    """Without out, the result is a new header; without alternative WCS, only the primary"""
    hdr_spec = fits.Header()
    hdr_spec["CTYPE1"] = "AWAV"
    hdr_spec["CRPIX1"] = 1
    hdr_spec["CRVAL1"] = 362.0
    hdr_spec["CDELT1"] = 1.86
    hdr_spec["CUNIT1"] = "nm"
    # the axis of the fibers of the RSS
    hdr_spec["CTYPE2"] = ""
    hdr_spec["CRPIX2"] = 0
    hdr_spec["CRVAL2"] = 0
    hdr_spec["CDELT2"] = 1
    hdr_sky = create_sky_header2()
    out = merge_wcs(hdr_sky, hdr_spec)

    assert out is not hdr_spec
    assert "CTYPE3" not in hdr_spec
    wcs0 = astropy.wcs.WCS(header=out, fix=False)
    assert list(wcs0.wcs.ctype) == ["RA---TAN", "DEC--TAN", "AWAV"]
    assert "CTYPE1B" not in out


@pytest.fixture
def rss_file(tmp_path, fiber_header):
    rss = make_rss(np.ones((623, 3)), fiber_header)
    rss[0].header["IPA"] = 30.0
    path = tmp_path / "rss.fits"
    rss.writeto(path)
    return path


@pytest.mark.parametrize("options", [[], ["-m", "linear"], ["-d"], ["--fix-missing"]])
def test_main(tmp_path, rss_file, options):
    outfile = tmp_path / "cube.fits"
    main([str(rss_file), "-o", str(outfile), "-p", "1.0"] + options)
    with fits.open(outfile) as cube:
        assert cube[0].data.shape[0] == 3
        assert cube[0].header["CDELT2"] == pytest.approx(1.0 / 3600)


def test_main_pa(tmp_path, rss_file, fiber_header):
    """The rotation is computed from IPA, or taken from the header with --wcs-pa-from-header"""
    outfile = tmp_path / "cube.fits"
    main([str(rss_file), "-o", str(outfile), "-p", "1.0"])
    with fits.open(outfile) as cube:
        pa = math.radians(compute_pa_from_ipa(30.0))
        assert cube[0].header["PC1_1"] == pytest.approx(math.cos(pa))
        assert cube[0].header["PC2_1"] == pytest.approx(math.sin(pa))

    main([str(rss_file), "-o", str(outfile), "-p", "1.0", "--wcs-pa-from-header"])
    with fits.open(outfile) as cube:
        assert cube[0].header["PC1_1"] == pytest.approx(fiber_header["PC1_1"])
        assert cube[0].header["PC2_1"] == pytest.approx(fiber_header["PC2_1"])
