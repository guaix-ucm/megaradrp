import pytest

from megaradrp.testing.create_tracemap import create_test_tracemap, create_test_tracemap2
from megaradrp.testing.create_wavecalib import create_test_wavecalib


@pytest.fixture
def wavecalib_data_state():
    return create_test_wavecalib()


@pytest.fixture
def tracemap_data_state():
    return create_test_tracemap()


@pytest.fixture
def tracemap_data():
    return create_test_tracemap2()