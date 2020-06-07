from skyfield.api import Angle
from skyfield.trigonometry import position_angle_of

def test_position_angle():
    a = Angle(degrees=0), Angle(degrees=0)
    b = Angle(degrees=1), Angle(degrees=1)
    assert str(position_angle_of(a, b)) == '315deg 00\' 15.7"'
