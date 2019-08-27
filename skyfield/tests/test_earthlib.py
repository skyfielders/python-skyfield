from numpy import array, sqrt
from skyfield.earthlib import AU_M, ERAD, reverse_terra, tau

def test_reverse_terra_with_zero_iterations():
    # With zero iterations, should return "geocentric" rather than
    # "geodetic" (="correct") longitude and latitude.
    lat, lon, elevation = reverse_terra(array([1, 0, 1]), 0, iterations=0)
    assert abs(lat - tau / 8) < 1e-16
    assert lon == 0.0
    assert abs(elevation - (AU_M * sqrt(2) - ERAD)) < 1e-16
