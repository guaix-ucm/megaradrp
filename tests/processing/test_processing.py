import os

import astropy.io.fits as fits
import pytest

from megaradrp.processing.trimover import trimOut
from megaradrp.processing.trimover import apextract_weights

DETCONF = {
    "trim1": [[0, 2056], [50, 4146]],
    "trim2": [[2156, 4212], [50, 4146]],
    "bng": [1, 1],
}


@pytest.fixture
def flat(simulated_flats):
    """A simulated flat, and the name of its file"""
    filename = simulated_flats.low[0]
    return fits.getdata(filename), str(filename)


@pytest.mark.parametrize("direction", ["normal", "mirror"])
def test_trim_and_o(flat, tmp_path, direction):
    data, filename = flat
    trimOut(filename, DETCONF, out=str(tmp_path / "result.fits"), direction=direction)
    with fits.open(tmp_path / "result.fits") as hdul:
        assert hdul[0].shape[0] + 100 == data.shape[0]
        assert hdul[0].shape[1] + 100 == data.shape[1]


def test_trim_and_o_fail(flat, tmp_path):
    _, filename = flat
    direction = "fails"
    with pytest.raises(ValueError) as excinfo:
        trimOut(filename, DETCONF, out=str(tmp_path / "result.fits"), direction=direction)
    assert excinfo.value.args[0] == f"{direction} must be either 'normal' or 'mirror'"


def test_trim_and_o_fail2(flat, tmp_path):
    _, filename = flat
    bins = "fail"
    detconf = dict(DETCONF, bng=bins)
    with pytest.raises(ValueError) as excinfo:
        trimOut(filename, detconf, out=str(tmp_path / "result.fits"))
    assert excinfo.value.args[0] == f"{bins} must be one if '11', '12', '21, '22'"


@pytest.mark.skipif(
    not os.path.exists("master_weights_LCB_10img_1exp.tar"),
    reason="no way of currently testing this without tar file",
)
def test_apextract_weights():
    import tarfile

    file_name = "master_weights_LCB_10img_1exp.tar"

    data = fits.getdata("fiberflat_frame.fits")
    rss = apextract_weights(data, tarfile.open(file_name, "r"))
    hdu_rss = fits.PrimaryHDU(rss)
    final = fits.HDUList([hdu_rss])
    final.writeto("rss.fits", overwrite=True)
    assert True
