import pytest
import numina.instrument.generic
from megaradrp.loader import load_drp
from megaradrp.testing.create_image import create_simple_frame


@pytest.mark.parametrize(
    "conf, uuix",
    [
        ["default", "9f3fecbf-376c-47ae-a188-313b7d829104"],
        [
            "9f3fecbf-376c-47ae-a188-313b7d829104",
            "9f3fecbf-376c-47ae-a188-313b7d829104",
        ],
        [
            "ca3558e3-e50d-4bbc-86bd-da50a0998a48",
            "ca3558e3-e50d-4bbc-86bd-da50a0998a48",
        ],
        [
            "9a86b2b2-3f7d-48ec-8f4f-3780ec967c90",
            "9a86b2b2-3f7d-48ec-8f4f-3780ec967c90",
        ],
        [
            "66f2283e-3049-4d4b-8ef1-14d62fcb611d",
            "66f2283e-3049-4d4b-8ef1-14d62fcb611d",
        ],
        [
            "4fd05b24-2ed9-457b-b563-a3c618bb1d4c",
            "4fd05b24-2ed9-457b-b563-a3c618bb1d4c",
        ],
    ],
)
def test_loader1(conf, uuix):
    import numina.core
    from numina.instrument.assembly import assembly_instrument

    obs = numina.core.ObservationResult(instrument="MEGARA")
    obs.frames.append(create_simple_frame(insconf=uuix))
    drpm = load_drp()
    obs.configuration = conf

    key, date_obs, keyname = drpm.select_profile(obs)
    ins = assembly_instrument(drpm.configurations, key, date_obs, by_key=keyname)
    assert isinstance(ins, numina.instrument.generic.InstrumentGeneric)
    assert str(ins.origin.uuid) == uuix
