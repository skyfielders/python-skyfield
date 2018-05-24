"""Tests of whether units behave."""

from assay import assert_raises
from numpy import array
from skyfield import units

try:
    from astropy import units as u
except ImportError:
    u = None

def needs_astropy(test):
    """Skip `test` if AstroPy is not available."""
    return None if (u is None) else test

def test_velocity_input_units():
    v1 = units.Velocity(au_per_d=2.0)
    v2 = units.Velocity(km_per_s=3462.9137)
    assert abs(v1.au_per_d - v2.au_per_d) < 1e-7

def test_stringifying_vector_distance():
    a = array([1.23, 4.56])
    s = str(units.Distance(au=a))
    if '[1' in s:
        # Python 3.5, says Travis CI.  No idea.
        assert s == '[1.23 4.56] au'
    else:
        # Every other version of Python.
        assert s == '[ 1.23  4.56] au'

def test_iterating_over_raw_measurement():
    distance = units.Distance(au=1.234)
    with assert_raises(units.UnpackingError):
        x, y, z = distance

def test_iterating_over_raw_velocity():
    velocity = units.Velocity(au_per_d=1.234)
    with assert_raises(units.UnpackingError):
        x, y, z = velocity

def test_converting_from_km_to_m():
    distance = units.Distance(km=1.234)
    assert abs(distance.m - 1234.0) < 1e-15

def test_converting_from_m_to_km():
    distance = units.Distance(m=1234.0)
    assert abs(distance.km - 1.234) < 1e-15

@needs_astropy
def test_converting_distance_with_astropy():
    distance = units.Distance(au=1.234)
    value1 = distance.km
    value2 = distance.to(u.km)
    epsilon = 0.02         # definitions of AU seem to disagree slightly
    assert abs(value1 - value2.value) < epsilon

@needs_astropy
def test_converting_velocity_with_astropy():
    velocity = units.Velocity(au_per_d=1.234)
    value1 = velocity.km_per_s
    value2 = velocity.to(u.km / u.s)
    epsilon = 1e-6
    assert abs(value1 - value2.value) < epsilon
