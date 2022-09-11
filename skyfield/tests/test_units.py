"""Tests of whether units behave."""

from assay import assert_raises
from numpy import array, nan
from skyfield.units import (
    Angle, Distance, Velocity, UnpackingError, WrongUnitError,
)

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
    assert str(Angle(degrees=tenth * -0.25)) == '00deg 00\' 00.0"'

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
    assert str(Angle(degrees=array(91), signed=True)) == '''+91deg 00' 00.0"'''
    assert str(Angle(hours=array(12))) == '''12h 00m 00.00s'''

def test_angle_array_strs():
    h = Angle(hours=array([0.5, nan, -13]))
    d = Angle(degrees=h._degrees)

    assert str(h) == '3 values from 00h 30m 00.00s to -13h 00m 00.00s'
    assert str(d) == '''3 values from 07deg 30' 00.0" to -195deg 00' 00.0"'''

    with assert_raises(WrongUnitError):
        h.dstr()
        d.hstr()

    assert h.hstr() == d.hstr(warn=False) == [
        '00h 30m 00.00s',
        'nan',
        '-13h 00m 00.00s',
    ]
    assert d.dstr() == h.dstr(warn=False) == [
        '07deg 30\' 00.0"',
        'nan',
        '-195deg 00\' 00.0"',
    ]

    empty = Angle(radians=[])
    assert str(empty) == 'Angle []'
    assert empty.hstr(warn=False) == []
    assert empty.dstr() == []

    assert h.hstr(format='{0} {1} {2} {3} {4} {5}', places=6) == [
        ' 0 30 0 0 6', 'nan', '- 13 0 0 0 6']
    assert d.dstr(format='{0} {1} {2} {3} {4} {5}', places=6) == [
        ' 7 30 0 0 6', 'nan', '- 195 0 0 0 6']

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

def test_distance_input_units():
    for d in (
        Distance(1.0),
        Distance(au=1.0),       # deprecated
        Distance(m=149597870700),
        Distance(km=149597870.700),
        Distance.au(1.0),       # modern
        Distance.m(149597870700),
        Distance.km(149597870.700),
    ):
        assert abs(d.au - 1.0) <= 0

def test_velocity_input_units():
    v1 = Velocity(au_per_d=2.0)
    v2 = Velocity(km_per_s=3462.9137)
    assert abs(v1.au_per_d - v2.au_per_d) < 1e-7

    v1 = Velocity.au_per_d(2.0)
    v2 = Velocity.km_per_s(3462.9137)
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

def test_helpful_exceptions():
    distance = Distance(1.234)
    expect = '''\
to use this Distance, ask for its value in a particular unit:

    distance.au
    distance.km
    distance.m'''

    with assert_raises(UnpackingError) as a:
        x, y, z = distance
    assert str(a.exception) == expect

    with assert_raises(UnpackingError) as a:
        distance[0]
    assert str(a.exception) == expect

    velocity = Velocity(1.234)
    expect = '''\
to use this Velocity, ask for its value in a particular unit:

    velocity.au_per_d
    velocity.km_per_s
    velocity.m_per_s'''

    with assert_raises(UnpackingError) as a:
        x, y, z = velocity
    assert str(a.exception) == expect

    with assert_raises(UnpackingError) as a:
        velocity[0]
    assert str(a.exception) == expect

    angle = Angle(radians=1.234)
    expect = '''\
to use this Angle, ask for its value in a particular unit:

    angle.degrees
    angle.hours
    angle.radians'''

    with assert_raises(UnpackingError) as a:
        x, y, z = angle
    assert str(a.exception) == expect

    with assert_raises(UnpackingError) as a:
        angle[0]
    assert str(a.exception) == expect

def test_constructors_accept_plain_lists():
    Distance(au=[1,2,3])
    Distance(km=[1,2,3])
    Distance(m=[1,2,3])
    Velocity(au_per_d=[1,2,3])
    Velocity(km_per_s=[1,2,3])

def test_converting_from_km_to_m():
    distance = Distance(km=1.234)
    assert abs(distance.m - 1234.0) < 1e-15

def test_converting_from_m_to_km():
    distance = Distance(m=1234.0)
    assert abs(distance.km - 1.234) < 1e-15

def test_deprecated_method_from_au():
    distance = Distance.from_au(1.25)
    assert distance.au == 1.25

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
