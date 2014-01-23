"""Basic tests of the Skyfield API module and its contents."""

from skyfield import api

def test_whether_planets_have_radii():
    assert api.mercury.radius.km == 2440.0
    for planet in api.nine_planets:
        assert planet.radius.km > 0.0
