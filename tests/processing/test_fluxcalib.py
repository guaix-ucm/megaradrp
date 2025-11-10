import pytest

import numpy

from megaradrp.processing.fluxcalib import update_flux_limits, WAVLIM_KEYS, PIXLIM_KEYS
from megaradrp.testing.create_wcs import generate_wcs


@pytest.mark.parametrize("wcsl", [None, generate_wcs()])
@pytest.mark.parametrize("ref", [0, 1])
def test_update_flux_limits2(wcsl, ref):

    header = {}

    values = numpy.array([1, 1000, 2, 900, 100, 800])

    off = (ref + 1) % 2

    if wcsl is None:
        results = values + off
    else:
        lm = numpy.array([values, numpy.zeros_like(values)])
        wavelen_ = wcsl.all_pix2world(lm.T, ref)
        results = wavelen_[:, 0]

    fluxlimits = dict(zip(PIXLIM_KEYS, values))
    hdr = update_flux_limits(header, fluxlimits, wcs=wcsl, ref=ref)

    for key in WAVLIM_KEYS:
        assert key in hdr

    for key in PIXLIM_KEYS:
        assert key in hdr

    for key, val in zip(PIXLIM_KEYS, values):
        assert hdr[key][0] == val + off

    for key, res in zip(WAVLIM_KEYS, results):
        assert hdr[key][0] == res


def test_update_flux_invalid():

    header = {}

    values = [1, 1000, 2, 900, 100, 800]
    fluxlimits = dict(zip(PIXLIM_KEYS, values))

    ref = 2
    with pytest.raises(ValueError):
        update_flux_limits(header, fluxlimits, ref=ref)
