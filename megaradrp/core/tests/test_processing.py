import shutil
from tempfile import mkdtemp

import astropy.io.fits as fits
import numpy as np
import pytest

from megaradrp.core.processing import trim_and_o
from megaradrp.core.processing import apextract_weights
from megaradrp.simulation.detector import ReadParams, MegaraDetectorSat
from megaradrp.simulation.actions import simulate_flat


def generate_bias_file():
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    ron = 2.0
    gain = 1.0
    bias = 1000.0

    qe = 0.8 * np.ones(DSHAPE)
    qe[0:15, 0:170] = 0.0

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, qe=qe,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    return simulate_flat(detector, exposure=1.0, source=5000.0)

@pytest.mark.parametrize("direction", ['normal', 'mirror'])
def test_trim_and_o(direction):
    temporary_path = mkdtemp()
    fs = generate_bias_file()
    fits.writeto('%s/flat.fits' % (temporary_path), fs, clobber=True)
    trim_and_o('%s/flat.fits' % (temporary_path), out='%s/result.fits' % (temporary_path), direction=direction)
    with fits.open('%s/result.fits' % (temporary_path)) as hdul:
        assert hdul[0].shape[0] + 100 == fs.shape[0]
        assert hdul[0].shape[1] + 100 == fs.shape[1]

    shutil.rmtree(temporary_path)

def test_trim_and_o_fail():
    temporary_path = mkdtemp()
    fs = generate_bias_file()
    fits.writeto('%s/flat.fits' % (temporary_path), fs, clobber=True)

    direction = 'fails'

    with pytest.raises(ValueError) as excinfo:
        trim_and_o('%s/flat.fits' % (temporary_path), out='%s/result.fits' % (temporary_path), direction=direction)
    shutil.rmtree(temporary_path)
    assert excinfo.value.args[0] == "%s must be either 'normal' or 'mirror'" % direction

def test_trim_and_o_fail2():
    temporary_path = mkdtemp()
    fs = generate_bias_file()
    fits.writeto('%s/flat.fits' % (temporary_path), fs, clobber=True)

    bins = 'fail'

    with pytest.raises(ValueError) as excinfo:
        trim_and_o('%s/flat.fits' % (temporary_path), out='%s/result.fits' % (temporary_path), bins=bins)
    shutil.rmtree(temporary_path)
    assert excinfo.value.args[0] == "%s must be one if '11', '12', '21, '22'" % bins


@pytest.mark.skip(reason="no way of currently testing this without tar file")
def test_apextract_weights():
    import tarfile

    file_name='master_weights_LCB_10img_1exp.tar'

    data = fits.getdata('fiberflat_frame.fits')
    rss = apextract_weights(data, tarfile.open(file_name, 'r'))
    hdu_rss = fits.PrimaryHDU(rss)
    final = fits.HDUList([hdu_rss])
    final.writeto('rss.fits', clobber=True)
    assert True



if __name__ == "__main__":
    # test_trim_and_o()
    # test_trim_and_o_fail()
    # test_trim_and_o_fail2()
    test_apextract_weights()