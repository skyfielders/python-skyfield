
from skyfield import magnitudelib as m
from skyfield.tests.conventions import A

def test_earth_magnitude_function():
    mag = m._earth_magnitude(0.983331936476, 1.41317594650699, 8.7897)
    assert abs(-3.269 - mag) < 0.0005
    mag = m._earth_magnitude(0.983356079811, 0.26526856764726, 4.1369)
    assert abs(-6.909 - mag) < 0.0005
    mag = m._earth_magnitude(0.983356467727, 0.62933287342927, 175.6869)
    assert abs(1.122 - mag) < 0.0005

    args = [
        A[0.983331936476, 0.983356079811, 0.983356467727],
        A[1.41317594650699, 0.26526856764726, 0.62933287342927],
        A[8.7897, 4.1369, 175.6869],
    ]
    magnitudes = m._earth_magnitude(*args)
    expected = [-3.269, -6.909, 1.122]
    assert all(magnitudes - expected < 0.0005)

def test_jupiter_magnitude_function():
    mag = m._jupiter_magnitude(5.446231815414, 6.44985867459088, 0.2446)
    assert abs(-1.667 - mag) < 0.0005
    mag = m._jupiter_magnitude(4.957681473205, 3.95393078136013, 0.3431)
    assert abs(-2.934 - mag) < 0.0005
    mag = m._jupiter_magnitude(5.227587855371, 5.23501920009381, 147.0989)
    assert abs(0.790 - mag) < 0.0005

    args = [
        A[5.446231815414, 4.957681473205, 5.227587855371],
        A[6.44985867459088, 3.95393078136013, 5.23501920009381],
        A[0.2446, 0.3431, 147.0989],
    ]
    magnitudes = m._jupiter_magnitude(*args)
    expected = [-1.667, -2.934, 0.790]
    assert all(magnitudes - expected < 0.0005)

def test_mercury_magnitude_function():
    mag = m._mercury_magnitude(0.310295423552, 1.32182643625754, 1.1677)
    assert abs(-2.477 - mag) < 0.0005
    mag = m._mercury_magnitude(0.413629222334, 0.92644808718613, 90.1662)
    assert abs(0.181 - mag) < 0.0005
    mag = m._mercury_magnitude(0.448947624811, 0.56004973217883, 178.7284)
    assert abs(7.167 - mag) < 0.0005

    args = [
        A[0.310295423552, 0.413629222334, 0.448947624811],
        A[1.32182643625754, 0.92644808718613, 0.56004973217883],
        A[1.1677, 90.1662, 178.7284],
    ]
    magnitudes = m._mercury_magnitude(*args)
    expected = [-2.477, 0.181, 7.167]
    assert all(magnitudes - expected < 0.0005)

def test_uranus_magnitude_function():
    mag = m._uranus_magnitude(18.321003215845, 17.3229728525108, 0.0410, -20.29, -20.28)
    assert abs(5.381 - mag) < 0.0005
    mag = m._uranus_magnitude(20.096361095266, 21.0888470145276, 0.0568, 1.02, 0.97)
    assert abs(6.025 - mag) < 0.0005
    mag = m._uranus_magnitude(19.38003071775, 11.1884243801383, 161.7728, -71.16, 55.11)
    assert abs(8.318 - mag) < 0.0005

    args = [
        A[18.321003215845, 20.096361095266, 19.38003071775],
        A[17.3229728525108, 21.0888470145276, 11.1884243801383],
        A[0.0410, 0.0568, 161.7728],
        A[-20.29, 1.02, -71.16],
        A[-20.28, 0.97, 55.11],
    ]
    magnitudes = m._uranus_magnitude(*args)
    expected = [5.381, 6.025, 8.318]
    assert all(magnitudes - expected < 0.0005)

def test_venus_magnitude_function():
    mag = m._venus_magnitude(0.722722540169, 1.71607489554051, 1.3232)
    assert abs(-3.917 - mag) < 0.0005
    mag = m._venus_magnitude(0.721480714554, 0.37762511206278, 124.1348)
    assert abs(-4.916 - mag) < 0.0005
    mag = m._venus_magnitude(0.726166592736, 0.28889582420642, 179.1845)
    assert abs(-3.090 - mag) < 0.0005

    args = [
        A[0.722722540169, 0.721480714554, 0.726166592736],
        A[1.71607489554051, 0.37762511206278, 0.28889582420642],
        A[1.3232, 124.1348, 179.1845],
    ]
    magnitudes = m._venus_magnitude(*args)
    expected = [-3.917, -4.916, -3.090]
    assert all(magnitudes - expected < 0.0005)
