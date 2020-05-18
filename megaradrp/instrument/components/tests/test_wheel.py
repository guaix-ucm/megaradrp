import pytest

from ..wheel import VPHWheel


class Vph(object):
    def __init__(self, name):
        self.name = name


@pytest.fixture
def wheel_dev():
    wheel = VPHWheel(3)
    for idx in range(3):
        wheel.put_in_pos(Vph(idx), idx)
    return wheel


def test_wheel(wheel_dev):
    curr = wheel_dev.current()
    assert isinstance(curr, Vph)
    assert curr.name == 0
    assert wheel_dev.pos() == 0