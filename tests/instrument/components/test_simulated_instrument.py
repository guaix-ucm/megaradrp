"""Components of the simulated instrument"""

import types

import astropy.units as u
import numpy as np
import pytest

import megaradrp.instrument.components.instrument as ins
from megaradrp.instrument.components.cover import HemiCover
from megaradrp.instrument.components.lamps import BlackBodyLamp


def test_blackbody_lamp():
    """The flux per unit wavelength of a black body, by Planck's law"""
    import astropy.constants as const

    wl = [5000.0, 8000.0] * u.AA
    temp = 3000 * u.K
    flux = BlackBodyLamp("BB", temp, factor=2.0).flux(wl)
    expected = 2 * 2 * const.h * const.c**2 / wl**5 / (np.exp(const.h * const.c / (wl * const.k_B * temp)) - 1) / u.sr
    assert flux.unit.is_equivalent(u.erg / u.s / u.cm**2 / u.AA / u.sr)
    assert np.allclose(flux.to_value(expected.unit), expected.value)


def test_hemicover():
    cover = HemiCover("cover")
    opened = []
    cover.opened.connect(lambda: opened.append(True))
    cover.open()
    assert cover.position == 1
    assert opened == [True]


def instrument(nfibers=3):
    """A MegaraInstrument with fakes of the components used in the projection"""
    obj = ins.MegaraInstrument.__new__(ins.MegaraInstrument)
    one = types.SimpleNamespace(transmission=lambda wl: 1.0, qe_wl=lambda wl: 1.0)
    obj.fiberset = types.SimpleNamespace(nfibers=nfibers, sigma=1.5, transmission=one.transmission)
    obj.internal_optics = obj.vph = one
    obj.detector = types.SimpleNamespace(qe_wl=one.qe_wl)
    obj.pseudo_slit = "pseudo_slit"
    obj.focal_actuator = types.SimpleNamespace(internal_focus_factor=2.0)
    obj.focal_plane = types.SimpleNamespace(get_visible_fibers=lambda: [(1, None), (3, None)])
    obj.get_visible_fibers = lambda: {"fibid": np.array([1, 3]), "cover": np.array([1.0, 0.5])}
    return obj


@pytest.fixture
def projections(monkeypatch):
    """The arguments of the projections of the fibers in the detector"""
    calls = []

    def record(name):
        def projection(*args):
            calls.append((name, args))
            return name

        return projection

    monkeypatch.setattr(ins, "project_rss", record("project_rss"))
    monkeypatch.setattr(ins, "project_rss_w", record("project_rss_w"))
    return calls


def test_simulate_focal_plane(projections):
    obj = instrument()
    wl = np.array([5000.0, 6000.0])
    assert obj.simulate_focal_plane(wl, np.ones_like(wl)) == "project_rss"
    ((_, args),) = projections
    fibids, pseudo_slit, vph, detector, sigma, wl_in, spec_in = args
    assert list(fibids) == [1, 3]
    # the sigma of the fibers with the focus factor
    assert sigma == 3.0
    # the photons of the visible fibers, with their coverage
    assert np.allclose(spec_in[:, 0], [1.0, 0.0, 0.5])


def test_project_rss_w(projections):
    obj = instrument()
    assert obj.project_rss_w() == "project_rss_w"
    ((_, args),) = projections
    assert args == ([1, 3], obj.pseudo_slit, obj.vph, obj.detector, 1.5)
