import pytest
import astropy.io.fits as fits
import numina.instrument.generic
import numina.core
from megaradrp.loader import load_drp


def create_simple_frame():
    hdr = {'INSTRUME': 'MEGARA', 'INSCONF': 'ca3558e3-e50d-4bbc-86bd-da50a0998a48'}
    hdu = fits.PrimaryHDU(data=[[1]])
    for k, v in hdr.items():
        hdu.header[k] = v
    hdulist = fits.HDUList([hdu])
    frame = numina.core.DataFrame(frame=hdulist)
    return frame


@pytest.mark.parametrize("conf, uuix", [
    ['default', "ca3558e3-e50d-4bbc-86bd-da50a0998a48"],
    ["ca3558e3-e50d-4bbc-86bd-da50a0998a48", "ca3558e3-e50d-4bbc-86bd-da50a0998a48"],
    ["9a86b2b2-3f7d-48ec-8f4f-3780ec967c90", "9a86b2b2-3f7d-48ec-8f4f-3780ec967c90"],
    ["66f2283e-3049-4d4b-8ef1-14d62fcb611d", "66f2283e-3049-4d4b-8ef1-14d62fcb611d"],
    ["4fd05b24-2ed9-457b-b563-a3c618bb1d4c", "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"]
])
def test_loader1(conf, uuix):
    import numina.core
    from numina.instrument.assembly import assembly_instrument

    obs = numina.core.ObservationResult(instrument='MEGARA')
    obs.frames.append(create_simple_frame())
    drpm = load_drp()
    obs.configuration = conf

    key, date_obs, keyname = drpm.select_profile(obs)
    ins = assembly_instrument(drpm.configurations, key, date_obs, by_key=keyname)
    assert isinstance(ins, numina.instrument.generic.InstrumentGeneric)
    assert str(ins.origin.uuid) == uuix
