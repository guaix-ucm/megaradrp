import pytest

from ..cube import create_cube


def test_create_cube_raise():
    with pytest.raises(ValueError):
        create_cube(None, None, 3)