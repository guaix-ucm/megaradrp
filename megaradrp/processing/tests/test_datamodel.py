
import pytest
import astropy.io.fits as fits

from ..datamodel import MegaraDataModel, FibersConf


def create_empty_img(insmode):
    img = fits.HDUList([fits.PrimaryHDU()])
    img[0].header['insmode'] = insmode
    return img


def test_fiberconf():

    datamodel = MegaraDataModel()

    img = create_empty_img('LCB')

    conf = datamodel.get_fiberconf(img)

    assert isinstance(conf, FibersConf)
    # Default values from file
    assert conf.name == 'LCB'
    assert conf.conf_id == 1
    assert conf.nbundles == 89
    assert conf.nfibers ==  623


def test_fiberconf_empty_mos():

    datamodel = MegaraDataModel()

    img = create_empty_img('MOS')

    with pytest.raises(ValueError):
        datamodel.get_fiberconf(img)
