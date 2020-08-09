from skyfield import magnitudelib as m
from skyfield.api import load

def test_front_end_function():
    # Simply call the routine with each planet to discover any exceptions.
    ts = load.timescale()
    t = ts.utc(2020, 7, 31)
    eph = load('de421.bsp')
    for name in ('mercury', 'venus', 'earth',
                 'jupiter barycenter', 'uranus barycenter'):
        astrometric = eph['sun'].at(t).observe(eph[name])
        m.planetary_magnitude(astrometric)

def test_mercury_magnitude_function():
    assert abs(-2.477 - m._mercury_magnitude(0.310295423552, 1.32182643625754, 1.1677)) < 0.0005
    assert abs(0.181 - m._mercury_magnitude(0.413629222334, 0.92644808718613, 90.1662)) < 0.0005
    assert abs(7.167 - m._mercury_magnitude(0.448947624811, 0.56004973217883, 178.7284)) < 0.0005

def test_venus_magnitude_function():
    assert abs(-3.917 - m._venus_magnitude(0.722722540169, 1.71607489554051, 1.3232)) < 0.0005
    assert abs(-4.916 - m._venus_magnitude(0.721480714554, 0.37762511206278, 124.1348)) < 0.0005
    assert abs(-3.090 - m._venus_magnitude(0.726166592736, 0.28889582420642, 179.1845)) < 0.0005

def test_earth_magnitude_function():
    assert abs(-3.269 - m._earth_magnitude(0.983331936476, 1.41317594650699, 8.7897)) < 0.0005
    assert abs(-6.909 - m._earth_magnitude(0.983356079811, 0.26526856764726, 4.1369)) < 0.0005
    assert abs(1.122 - m._earth_magnitude(0.983356467727, 0.62933287342927, 175.6869)) < 0.0005

def test_mars_magnitude_function():
    pass

def test_jupiter_magnitude_function():
    assert abs(-1.667 - m._jupiter_magnitude(5.446231815414, 6.44985867459088, 0.2446)) < 0.0005
    assert abs(-2.934 - m._jupiter_magnitude(4.957681473205, 3.95393078136013, 0.3431)) < 0.0005
    assert abs(0.790 - m._jupiter_magnitude(5.227587855371, 5.23501920009381, 147.0989)) < 0.0005

def test_saturn_magnitude_function():
    pass

def test_uranus_magnitude_function():
    assert abs(5.381 - m._uranus_magnitude(18.321003215845, 17.3229728525108, 0.0410, -20.29, -20.28)) < 0.0005
    assert abs(6.025 - m._uranus_magnitude(20.096361095266, 21.0888470145276, 0.0568, 1.02, 0.97)) < 0.0005
    assert abs(8.318 - m._uranus_magnitude(19.38003071775, 11.1884243801383, 161.7728, -71.16, 55.11)) < 0.0005

def test_neptune_magnitude_function():
    pass
