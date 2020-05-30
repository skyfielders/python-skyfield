"""Tests of whether units behave."""

from assay import assert_raises
from numpy import array
from skyfield.units import Angle, Distance, Velocity, UnpackingError

try:
    from astropy import units as u
except ImportError:
    u = None

def needs_astropy(test):
    """Skip `test` if AstroPy is not available."""
    return None if (u is None) else test

def test_degree_rounding():
    tenth = 0.1 / 60.0 / 60.0  # of an arcsecond

    assert str(Angle(degrees=tenth * -600.75)) == '-00deg 01\' 00.1"'
    assert str(Angle(degrees=tenth * -600.25)) == '-00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * -599.75)) == '-00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * -599.25)) == '-00deg 00\' 59.9"'

    assert str(Angle(degrees=tenth * -1.75)) == '-00deg 00\' 00.2"'
    assert str(Angle(degrees=tenth * -1.25)) == '-00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * -0.75)) == '-00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * -0.25)) == '-00deg 00\' 00.0"'

    assert str(Angle(degrees=0.0)) == '00deg 00\' 00.0"'

    assert str(Angle(degrees=tenth * 0.25)) == '00deg 00\' 00.0"'
    assert str(Angle(degrees=tenth * 0.75)) == '00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * 1.25)) == '00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * 1.75)) == '00deg 00\' 00.2"'

    assert str(Angle(degrees=tenth * 599.25)) == '00deg 00\' 59.9"'
    assert str(Angle(degrees=tenth * 599.75)) == '00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * 600.25)) == '00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * 600.75)) == '00deg 01\' 00.1"'

def test_angle_scalar_strs():
    assert str(Angle(degrees=array(91))) == '''91deg 00' 00.0"'''
    assert str(Angle(hours=array(12))) == '''12h 00m 00.00s'''

def test_angle_array_strs():
    assert str(Angle(degrees=array([90, 91, 92]))) == (
        '''3 values from 90deg 00' 00.0" to 92deg 00' 00.0"'''
        )
    assert str(Angle(hours=array([11, 12, 13]))) == (
        '''3 values from 11h 00m 00.00s to 13h 00m 00.00s'''
        )

def test_angle_sexagesimal_args():
    assert str(Angle(degrees=(90,))) == '''90deg 00' 00.0"'''
    assert str(Angle(hours=(12,))) == '''12h 00m 00.00s'''

    assert str(Angle(degrees=(90, 15))) == '''90deg 15' 00.0"'''
    assert str(Angle(hours=(12, 30))) == '''12h 30m 00.00s'''

    assert str(Angle(degrees=(90, 15, 30))) == '''90deg 15' 30.0"'''
    assert str(Angle(hours=(12, 30, 15))) == '''12h 30m 15.00s'''

def test_arcminutes_and_arcseconds_and_mas():
    angle = Angle(degrees=1.0)
    assert angle.arcminutes() == 60
    assert angle.arcseconds() == 60 * 60
    assert angle.mas() == 60 * 60 * 1000

def test_velocity_input_units():
    v1 = Velocity(au_per_d=2.0)
    v2 = Velocity(km_per_s=3462.9137)
    assert abs(v1.au_per_d - v2.au_per_d) < 1e-7

def test_stringifying_vector_distance():
    a = array([1.23, 4.56])
    s = str(Distance(au=a))
    if '[1' in s:
        # Python 3.5, says Travis CI.  No idea.
        assert s == '[1.23 4.56] au'
    else:
        # Every other version of Python.
        assert s == '[ 1.23  4.56] au'

def test_iterating_over_raw_measurement():
    distance = Distance(au=1.234)
    with assert_raises(UnpackingError) as a:
        x, y, z = distance
    assert str(a.exception) == '''\
cannot directly unpack a Distance into several values

To unpack a Distance into three components, you need to ask for its
value in specific units through an attribute or method:

    x, y, z = distance.au
    x, y, z = distance.km
    x, y, z = distance.to(astropy_unit)
'''

def test_iterating_over_raw_velocity():
    velocity = Velocity(au_per_d=1.234)
    with assert_raises(UnpackingError) as a:
        x, y, z = velocity
    assert str(a.exception) == '''\
cannot directly unpack a Velocity into several values

To unpack a Velocity into three components, you need to ask for its
value in specific units through an attribute or method:

    xdot, ydot, zdot = velocity.au_per_d
    xdot, ydot, zdot = velocity.km_per_s
    xdot, ydot, zdot = velocity.to(astropy_unit)
'''

def test_iterating_over_raw_angle():
    angle = Angle(degrees=4.5)
    with assert_raises(ValueError) as a:
        iter(angle)
    assert str(a.exception) == '''choose a specific Angle unit to iterate over

Instead of iterating over this Angle object, try iterating over one of
its unit-specific arrays like .degrees, .hours, or .radians, or else over
the output of one of its methods like .hstr(), .dstr(), .arcminutes(),
.arcseconds(), or .mas().  For all of the possibilities see:
https://rhodesmill.org/skyfield/api-units.html#skyfield.units.Angle'''

def test_converting_from_km_to_m():
    distance = Distance(km=1.234)
    assert abs(distance.m - 1234.0) < 1e-15

def test_converting_from_m_to_km():
    distance = Distance(m=1234.0)
    assert abs(distance.km - 1.234) < 1e-15

@needs_astropy
def test_converting_distance_with_astropy():
    distance = Distance(au=1.234)
    value1 = distance.km
    value2 = distance.to(u.km)
    epsilon = 0.02         # definitions of AU seem to disagree slightly
    assert abs(value1 - value2.value) < epsilon

@needs_astropy
def test_converting_velocity_with_astropy():
    velocity = Velocity(au_per_d=1.234)
    value1 = velocity.km_per_s
    value2 = velocity.to(u.km / u.s)
    epsilon = 1e-6
    assert abs(value1 - value2.value) < epsilon
