import textwrap
from ..api import Topos, load, wgs84
from ..sgp4lib import EarthSatellite

lines = [
    'ISS (ZARYA)             ',
    '1 25544U 98067A   13330.58127943  .00000814  00000-0  21834-4 0  1064',
    '2 25544  51.6484  23.7537 0001246  74.1647  18.7420 15.50540527859894',
]

def dedent(s):
    return textwrap.dedent(s.rstrip())

def eph():
    yield load('de421.bsp')

def test_jpl_segment(eph):
    e = eph['mercury barycenter']
    expected = dedent("""\
        'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 1 MERCURY BARYCENTER
    """)
    assert str(e) == expected
    expected = dedent("""\
        <ChebyshevPosition 'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 1 MERCURY BARYCENTER>
    """)
    assert repr(e) == expected

def test_satellite_with_name():
    s = EarthSatellite(lines[1], lines[2], lines[0])
    expected = dedent("""\
        ISS (ZARYA) catalog #25544 epoch 2013-11-26 13:57:03 UTC
    """)
    assert str(s) == expected
    expected = dedent("""\
        <EarthSatellite ISS (ZARYA) catalog #25544 epoch 2013-11-26 13:57:03 UTC>
    """)
    assert repr(s) == expected

def test_satellite_without_name():
    s = EarthSatellite(lines[1], lines[2])
    expected = dedent("""\
        catalog #25544 epoch 2013-11-26 13:57:03 UTC
    """)
    assert str(s) == expected
    expected = dedent("""\
        <EarthSatellite catalog #25544 epoch 2013-11-26 13:57:03 UTC>
    """)
    assert repr(s) == expected

def test_geographic_position():
    t = wgs84.latlon(42.2, -88.1)
    expected = dedent("""\
        WGS84 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m
    """)
    assert str(t) == expected
    expected = dedent("""\
        <GeographicPosition WGS84 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m>
    """)
    assert repr(t) == expected

    t.target_name = 'Custom name'
    assert str(t) == 'Custom name'
    assert repr(t) == "<GeographicPosition Custom name>"

def test_topos():
    t = Topos(latitude_degrees=42.2, longitude_degrees=-88.1)
    expected = dedent("""\
        IERS2010 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m
    """)
    assert str(t) == expected
    expected = dedent("""\
        <Topos IERS2010 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m>
    """)
    assert repr(t) == expected

    t.target_name = 'Custom name'
    assert str(t) == 'Custom name'
    assert repr(t) == "<Topos Custom name>"

def test_jpl_vector_sum(eph):
    e = eph['earth']
    expected = dedent("""\
        Sum of 2 vectors:
         'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
         'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH
    """)
    assert str(e) == expected
    expected = dedent("""\
        <VectorSum of 2 vectors:
         'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
         'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH>
    """)
    assert repr(e) == expected


def test_geographic_position_and_earth_satellite_vector_sum(eph):
    s = EarthSatellite(lines[1], lines[2])
    t = wgs84.latlon(42.2, -88.1)
    v = s - t
    expected = dedent("""\
        Sum of 2 vectors:
         Reversed Geodetic WGS84 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m -> 399 EARTH
         EarthSatellite 399 EARTH -> catalog #25544 epoch 2013-11-26 13:57:03 UTC
    """)
    assert str(v) == expected
    expected = dedent("""\
        <VectorSum of 2 vectors:
         Reversed Geodetic WGS84 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m -> 399 EARTH
         EarthSatellite 399 EARTH -> catalog #25544 epoch 2013-11-26 13:57:03 UTC>
    """)
    assert repr(v) == expected

def test_topos_and_earth_satellite_vector_sum(eph):
    s = EarthSatellite(lines[1], lines[2])
    t = Topos(latitude_degrees=42.2, longitude_degrees=-88.1)
    v = s - t
    expected = dedent("""\
        Sum of 2 vectors:
         Reversed Geodetic IERS2010 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m -> 399 EARTH
         EarthSatellite 399 EARTH -> catalog #25544 epoch 2013-11-26 13:57:03 UTC
    """)
    assert str(v) == expected
    expected = dedent("""\
        <VectorSum of 2 vectors:
         Reversed Geodetic IERS2010 latitude +42.2000 N longitude -88.1000 E elevation 0.0 m -> 399 EARTH
         EarthSatellite 399 EARTH -> catalog #25544 epoch 2013-11-26 13:57:03 UTC>
    """)
    assert repr(v) == expected
