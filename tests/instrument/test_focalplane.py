import pytest

import megaradrp.datamodel as dm


@pytest.fixture(scope="module")
def focalplane(request):
    name = request.param
    fp = dm.get_fiberconf_default(name)
    return fp


@pytest.mark.parametrize("focalplane", ["MOS"], indirect=["focalplane"])
@pytest.mark.parametrize(
    "fibid, around",
    [
        [372, [373, 374, 375]],
        [8, [9, 10, 11]],
        [623, [620, 621, 622]],
        [26, [23, 25, 28]],  # sky fiber
        [25, [22, 23, 24, 26, 27, 28]],  # sky fiber
        [487, [484, 485, 486, 488, 489, 490]],  # sky fiber
        [490, [487, 488, 489]],
    ],
)
def test_nearby_mos(focalplane, fibid, around):
    # Find the bundle the fiber belongs to
    result = focalplane.nearby_fibers(fibid)
    assert sorted(result) == sorted(around)


@pytest.mark.parametrize("focalplane", ["LCB"], indirect=["focalplane"])
@pytest.mark.parametrize(
    "fibid, around",
    [
        [372, [250, 376, 252, 388, 374, 238]],
        [8, [12, 10]],
        [623, [184, 185, 525, 529, 620, 622]],
        [26, [23, 25, 28]],  # sky fiber
        [25, [22, 23, 24, 26, 27, 28]],  # sky fiber
        [487, [484, 485, 486, 488, 489, 490]],  # sky fiber
        [490, [487, 488, 489]],
    ],
)
def test_nearby_lcb(focalplane, fibid, around):
    # Find the bundle the fiber belongs to
    result = focalplane.nearby_fibers(fibid)
    assert sorted(result) == sorted(around)


@pytest.mark.parametrize("focalplane", ["LCB"], indirect=["focalplane"])
@pytest.mark.parametrize("fibid", [700, -1, 624])
def test_nearby_exception_lcb(focalplane, fibid):
    with pytest.raises(ValueError):
        focalplane.nearby_fibers(fibid)


@pytest.mark.parametrize("focalplane", ["MOS"], indirect=["focalplane"])
@pytest.mark.parametrize("fibid", [700, -1, 645])
def test_nearby_exception_mos(focalplane, fibid):
    with pytest.raises(ValueError):
        focalplane.nearby_fibers(fibid)


@pytest.mark.parametrize("funit, warns", [(None, False), ("arcsec", True)])
def test_lcb_units(caplog, capsys, funit, warns):
    """A warning if the units of the fibers do not match the position of the reference fiber"""
    import logging

    from megaradrp.datamodel import create_default_fiber_header
    from megaradrp.instrument.focalplane import FocalPlaneConf

    hdr = create_default_fiber_header("LCB")
    if funit is not None:
        hdr["FUNIT"] = funit
    caplog.set_level(logging.WARNING)
    conf = FocalPlaneConf.from_header(hdr)
    assert conf.name == "LCB"
    assert ("is in mm, but the units of the fibers are arcsec" in caplog.text) == warns
    assert capsys.readouterr().out == ""
