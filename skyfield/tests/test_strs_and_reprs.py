import textwrap
from ..api import Topos, load
from ..sgp4lib import EarthSatellite

def dedent(s):
    return textwrap.dedent(s.rstrip())

def eph():
    yield load('de421.bsp')

def test_jpl_segment(eph):
    e = eph['mercury barycenter']
    expected = dedent("""\
        Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 1 MERCURY BARYCENTER
    """)
    assert str(e) == expected
    expected = dedent("""\
        <Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 1 MERCURY BARYCENTER>
    """)
    assert repr(e) == expected

def test_satellite_with_name(eph):
    lines = [
        'ISS (ZARYA)             ',
        '1 25544U 98067A   13330.58127943  .00000814  00000-0  21834-4 0  1064',
        '2 25544  51.6484  23.7537 0001246  74.1647  18.7420 15.50540527859894',
    ]
    s = EarthSatellite(lines[1], lines[2], lines[0])
    expected = dedent("""\
        EarthSatellite 'ISS (ZARYA)' number=25544 epoch=2013-11-26T13:57:03Z
    """)
    assert str(s) == expected
    expected = dedent("""\
        <EarthSatellite 'ISS (ZARYA)' number=25544 epoch=2013-11-26T13:57:03Z>
    """)
    assert repr(s) == expected

def test_satellite_without_name(eph):
    lines = [
        '1 25544U 98067A   13330.58127943  .00000814  00000-0  21834-4 0  1064',
        '2 25544  51.6484  23.7537 0001246  74.1647  18.7420 15.50540527859894',
    ]
    s = EarthSatellite(lines[0], lines[1])
    expected = dedent("""\
        EarthSatellite number=25544 epoch=2013-11-26T13:57:03Z
    """)
    assert str(s) == expected
    expected = dedent("""\
        <EarthSatellite number=25544 epoch=2013-11-26T13:57:03Z>
    """)
    assert repr(s) == expected

def test_topos(eph):
    t = Topos(latitude_degrees=42.2, longitude_degrees=-88.1)
    expected = dedent("""\
        Topos 42deg 12' 00.0" N -88deg 06' 00.0" E
    """)
    assert str(t) == expected
    expected = dedent("""\
        <Topos 42deg 12' 00.0" N -88deg 06' 00.0" E>
    """)
    assert repr(t) == expected

def test_vector_sum(eph):
    e = eph['earth']
    expected = dedent("""\
        Sum of 2 vectors:
         + Segment 'de421.bsp' 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
         + Segment 'de421.bsp' 3 EARTH BARYCENTER -> 399 EARTH
    """)
    assert str(e) == expected
    expected = dedent("""\
        <VectorSum of 2 vectors 0 SOLAR SYSTEM BARYCENTER -> 399 EARTH>
    """)
    assert repr(e) == expected
