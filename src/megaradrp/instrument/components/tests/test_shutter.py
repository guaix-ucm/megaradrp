import pytest

from ..shutter import MegaraShutter

@pytest.fixture
def shutter_dev():
    wheel = MegaraShutter()
    return wheel


def test_shutter1(shutter_dev):
    curr = shutter_dev.current()
    assert shutter_dev.pos() == 1
    assert curr.name == 'OPEN'


def test_shutter_open(shutter_dev):
    shutter_dev.open()
    curr = shutter_dev.current()
    assert shutter_dev.pos() == 1
    assert curr.name == 'OPEN'


def test_shutter_stop(shutter_dev):
    shutter_dev.close()
    curr = shutter_dev.current()
    assert shutter_dev.pos() == 0
    assert curr.name == 'STOP'
