from skyfield.api import Angle, Topos, load, load_file
from skyfield.trigonometry import position_angle_of

def test_position_angle():
    a = Angle(degrees=0), Angle(degrees=0)
    b = Angle(degrees=1), Angle(degrees=1)
    assert str(position_angle_of(a, b)) == '315deg 00\' 15.7"'

def test_position_angle_against_nasa_horizons():
    ts = load.timescale(builtin=True)
    t = ts.utc(2053, 10, 9)

    eph = load_file('./skyfield/tests/data/jup310-2053-10-08.bsp')
    boston = eph['earth'] + Topos(longitude_degrees=(-71, 3, 24.8),
                                  latitude_degrees=(42, 21, 24.1))

    b = boston.at(t)
    j = b.observe(eph['jupiter'])#.apparent()
    i = b.observe(eph['io'])#.apparent()

    a = position_angle_of(j.radec(epoch='date'), i.radec(epoch='date'))

    assert abs(a.degrees - 293.671) < 0.002
