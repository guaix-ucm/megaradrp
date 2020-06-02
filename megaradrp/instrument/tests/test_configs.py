import pytest
import numina.instrument.generic
from megaradrp.loader import load_drp


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
    drpm = load_drp()
    obs.configuration = conf

    key, date_obs, keyname = drpm.select_profile(obs)
    ins = assembly_instrument(drpm.configurations, key, date_obs, by_key=keyname)
    assert isinstance(ins, numina.instrument.generic.InstrumentGeneric)
    assert str(ins.origin.uuid) == uuix


@pytest.mark.parametrize("uuix, uuixc", [
    ["ca3558e3-e50d-4bbc-86bd-da50a0998a48",
     ['2e02e135-2325-47c9-9975-466b445b0b8b', '78f6d437-85d6-4d7a-bdb7-1e043368b442']
     ],
    ["9a86b2b2-3f7d-48ec-8f4f-3780ec967c90",
     ['2e02e135-2325-47c9-9975-466b445b0b8b', '97d48545-3258-473e-9311-00e390999d52']],
    ["66f2283e-3049-4d4b-8ef1-14d62fcb611d",
     ['2e02e135-2325-47c9-9975-466b445b0b8b', '715c7ced-f989-42a6-bb74-a5b5cda0495c']],
    ["4fd05b24-2ed9-457b-b563-a3c618bb1d4c",
     ['2e02e135-2325-47c9-9975-466b445b0b8b', '78f6d437-85d6-4d7a-bdb7-1e043368b442']]
])
def test_subconf(uuix, uuixc):
    import numina.core
    from numina.instrument.assembly import assembly_instrument

    obs = numina.core.ObservationResult(instrument='MEGARA')
    drpm = load_drp()
    obs.configuration = uuix

    key, date_obs, keyname = drpm.select_profile(obs)
    ins = assembly_instrument(drpm.configurations, key, date_obs, by_key=keyname)
    for c,uu in zip(ins.children.values(), uuixc):
        assert str(c.origin.uuid) == uu



def represent_comp(comp):
    print('name', comp.name)
    print('parent', comp.parent)
    print('properties', comp.properties)
    for c,k in comp.children.items():
        represent_comp(k)
    print('oname:', comp.origin.name)
    print('ouuid=', comp.origin.uuid)

    print('odate_start=', comp.origin.date_start)
    print('odate_end=', comp.origin.date_end)