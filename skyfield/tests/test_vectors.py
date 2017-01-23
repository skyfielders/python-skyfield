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
