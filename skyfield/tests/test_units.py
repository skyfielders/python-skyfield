"""Tests of whether units behave."""

from numpy import array
from skyfield import units

try:
    from astropy import units as u
except ImportError:
    u = None

def needs_astropy(test):
    """Skip `test` if AstroPy is not available."""
    return None if (u is None) else test

def test_stringifying_vector_distance():
    a = array([1.23, 4.56])
    assert str(units.Distance(au=a)) == '[ 1.23  4.56] au'

def test_converting_from_km_to_m():
    distance = units.Distance(km=1.234)
    assert abs(distance.m - 1234.0) < 1e-15

def test_converting_from_m_to_km():
    distance = units.Distance(m=1234.0)
    assert abs(distance.km - 1.234) < 1e-15
    
def test_indexing_distance():
    d = units.Distance(au=[1., 2., 3., 4., 5.])
    assert d.shape == (5,)
    d0 = d[0]
    assert d.au[0] == d0.au
    assert d.km[0] == d0.km
    assert d.m[0] == d0.m

def test_slicing_distance():
    d = units.Distance(au=[1., 2., 3., 4., 5.])
    assert d.shape == (5,)
    d24 = d[2:4]
    assert d24.shape == (2,)
    assert (d.au[2:4] == d24.au).all()
    assert (d.km[2:4] == d24.km).all()
    assert (d.m[2:4] == d24.m).all()
    
def test_indexing_velocity():
    v = units.Velocity(au_per_d=[1., 2., 3., 4., 5.])
    assert v.shape == (5,)
    v0 = v[0]
    assert v.au_per_d[0] == v0.au_per_d
    assert v.km_per_s[0] == v0.km_per_s

def test_slicing_velocity():
    v = units.Velocity(au_per_d=[1., 2., 3., 4., 5.])
    assert v.shape == (5,)
    v24 = v[2:4]
    assert v24.shape == (2,)
    assert (v.au_per_d[2:4] == v24.au_per_d).all()
    assert (v.km_per_s[2:4] == v24.km_per_s).all()
    
def test_indexing_angle():
    a = units.Angle(hours=array([1., 2., 3., 4., 5.]))
    assert a.shape == (5,)
    a0 = a[0]
    assert a.preference == a0.preference
    assert a.radians[0] == a0.radians
    assert a._degrees[0] == a0._degrees
    assert a._hours[0] == a0._hours

def test_slicing_angle():
    a = units.Angle(hours=array([1., 2., 3., 4., 5.]))
    assert a.shape == (5,)
    a24 = a[2:4]
    assert a24.shape == (2,)
    assert (a.radians[2:4] == a24.radians).all()
    assert (a._degrees[2:4] == a24._degrees).all()
    assert (a._hours[2:4] == a24._hours).all()

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
