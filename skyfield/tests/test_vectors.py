# Test the behavior of all combinations of vector.

from assay import assert_raises
from skyfield.api import Topos, load

def test_bad_addition():
    planets = load('de421.bsp')
    earth = planets['earth']
    mars = planets['mars']
    with assert_raises(ValueError, 'the center where the other vector starts'):
        earth + mars

def test_bad_subtraction():
    planets = load('de421.bsp')
    earth = planets['earth']
    usno = Topos('38.9215 N', '77.0669 W', elevation_m=92.0)
    with assert_raises(ValueError, 'if they both start at the same center'):
        earth - usno

def test_chebyshev_subtraction():
    planets = load('de421.bsp')
    v = planets['earth barycenter'] - planets['sun']

    assert str(v) == """\
Sum of 2 vectors:
 - Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 10 SUN
 + Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER"""

    assert repr(v) == "\
<VectorSum of 2 vectors 10 SUN -> 3 EARTH BARYCENTER>"

def test_vectors():
    ts = load.timescale()
    t = ts.tt(2017, 1, 23, 10, 44)

    planets = load('de421.bsp')
    earth = planets['earth']
    mars = planets['mars']

    v = earth

    assert str(v) == """\
Sum of 2 vectors:
 + Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 + Segment 'de421.bsp' 3 EARTH BARYCENTER -> 399 EARTH"""

    assert repr(v) == "\
<VectorSum of 2 vectors 0 SOLAR SYSTEM BARYCENTER -> 399 EARTH>"

    assert str(v.at(t)) == "\
<Barycentric position and velocity at date t center=0 target=399>"

    v = earth - mars

    assert str(v) == """\
Sum of 4 vectors:
 - Segment 'de421.bsp' 4 MARS BARYCENTER -> 499 MARS
 - Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 4 MARS BARYCENTER
 + Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 + Segment 'de421.bsp' 3 EARTH BARYCENTER -> 399 EARTH"""

    assert repr(v) == "\
<VectorSum of 4 vectors 499 MARS -> 399 EARTH>"

    assert str(v.at(t)) == "\
<Geometric position and velocity at date t center=499 target=399>"
