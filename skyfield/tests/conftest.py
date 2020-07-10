"""conftest.py to provide session wide fixtures"""


import pytest
from skyfield.api import load


@pytest.fixture(scope='session')
def ts():
    """Provide standard timescale for tests"""
    return load.timescale()
