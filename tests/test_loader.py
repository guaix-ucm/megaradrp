import numina.core.pipeline as ppl
from megaradrp.loader import load_drp


def test_loader():

    drpm = load_drp()
    assert isinstance(drpm, ppl.InstrumentDRP)
    assert drpm.name == "MEGARA"
