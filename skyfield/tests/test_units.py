"""Tests of whether units behave."""

import pytest
from skyfield import units

try:
    from astropy import units as u
except ImportError:
    u = None

needs_astropy = pytest.mark.skipif(u is None, reason='cannot import AstroPy')

def test_iterating_over_raw_measurement():
    distance = units.Distance(AU=1.234)
    with pytest.raises(units.UnpackingError):
        x, y, z = distance

def test_iterating_over_raw_velocity():
    velocity = units.Velocity(AU_per_d=1.234)
    with pytest.raises(units.UnpackingError):
        x, y, z = velocity

def test_converting_from_km_to_m():
    distance = units.Distance(km=1.234)
    assert abs(distance.m - 1234.0) < 1e-15

def test_converting_from_m_to_km():
    distance = units.Distance(m=1234.0)
    assert abs(distance.km - 1.234) < 1e-15

@needs_astropy
def test_converting_distance_with_astropy():
    distance = units.Distance(AU=1.234)
    value1 = distance.km
    value2 = distance.to(u.km)
    epsilon = 0.02         # definitions of AU seem to disagree slightly
    assert abs(value1 - value2.value) < epsilon

@needs_astropy
def test_converting_velocity_with_astropy():
    velocity = units.Velocity(AU_per_d=1.234)
    value1 = velocity.km_per_s
    value2 = velocity.to(u.km / u.s)
    epsilon = 1e-6
    assert abs(value1 - value2.value) < epsilon
