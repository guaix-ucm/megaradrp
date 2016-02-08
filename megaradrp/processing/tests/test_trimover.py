from megaradrp.processing.trimover import OverscanCorrector
from numina.core import DataFrame
import shutil
from tempfile import mkdtemp
import numpy as np
import astropy.io.fits as fits

def test_OverscanCorrector():

    oc = OverscanCorrector()
    assert oc.trim1 == (slice(0, 2056, None),slice(50, 4146, None))
    assert oc.pcol1 == (slice(0, 2056, None),slice(0, 50, None))
    assert oc.ocol1 == (slice(0, 2056, None),slice(4146, 4196, None))
    assert oc.orow1 == (slice(2056, 2106, None),slice(50, 4146, None))

    assert oc.trim2 == (slice(2156, 4212, None),slice(50, 4146, None))
    assert oc.pcol2 == (slice(2156, 4212, None),slice(4146, 4196, None))
    assert oc.ocol2 == (slice(2156, 4212, None),slice(0, 50, None))
    assert oc.orow2 == (slice(2106, 2156, None),slice(50, 4146, None))


if __name__ == "__main__":
    test_OverscanCorrector()
