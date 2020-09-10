# Test the behavior of all combinations of vector.

from assay import assert_raises
from skyfield.api import Topos, load
from skyfield.positionlib import Geocentric

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
 Reversed 'de421.bsp' segment 10 SUN -> 0 SOLAR SYSTEM BARYCENTER
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER"""

    assert repr(v) == """\
<VectorSum of 2 vectors:
 Reversed 'de421.bsp' segment 10 SUN -> 0 SOLAR SYSTEM BARYCENTER
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER>"""

def test_negation():
    ts = load.timescale()
    t = ts.utc(2020, 8, 30, 16, 5)
    usno = Topos('38.9215 N', '77.0669 W', elevation_m=92.0)
    neg = -usno
    p1 = usno.at(t)
    p2 = neg.at(t)
    assert (p1.position.au == - p2.position.au).all()
    assert (p1.velocity.au_per_d == - p2.velocity.au_per_d).all()

    # A second negation should return the unwrapped original.
    neg = -neg
    assert neg is usno

def test_vectors():
    ts = load.timescale()
    t = ts.tt(2017, 1, 23, 10, 44)

    planets = load('de421.bsp')
    earth = planets['earth']
    mars = planets['mars']

    v = earth

    assert str(v) == """\
Sum of 2 vectors:
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH"""

    assert repr(v) == """\
<VectorSum of 2 vectors:
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH>"""

    assert str(v.at(t)) == "\
<Barycentric BCRS position and velocity at date t center=0 target=399>"

    v = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    assert str(v) == """\
Sum of 3 vectors:
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH
 Topos 399 EARTH -> Earth latitude 38deg 55' 17.4" N longitude -77deg 04' 00.8" E elevation 92 m"""

    assert repr(v) == """\
<VectorSum of 3 vectors:
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH
 Topos 399 EARTH -> Earth latitude 38deg 55' 17.4" N longitude -77deg 04' 00.8" E elevation 92 m>"""

    assert str(v.at(t)) == """\
<Barycentric BCRS position and velocity at date t center=0 target=Earth latitude 38deg 55' 17.4" N longitude -77deg 04' 00.8" E elevation 92 m>"""

    v = earth - mars

    assert str(v) == """\
Sum of 4 vectors:
 Reversed 'de421.bsp' segment 499 MARS -> 4 MARS BARYCENTER
 Reversed 'de421.bsp' segment 4 MARS BARYCENTER -> 0 SOLAR SYSTEM BARYCENTER
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH"""

    assert repr(v) == """\
<VectorSum of 4 vectors:
 Reversed 'de421.bsp' segment 499 MARS -> 4 MARS BARYCENTER
 Reversed 'de421.bsp' segment 4 MARS BARYCENTER -> 0 SOLAR SYSTEM BARYCENTER
 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
 'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH>"""

    assert str(v.at(t)) == "\
<ICRF position and velocity at date t center=499 target=399>"

    geocentric = Geocentric([0,0,0])
    assert geocentric.center == 399
