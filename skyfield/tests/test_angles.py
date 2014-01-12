from skyfield.angles import Angle, HourAngle

def test_degree_rounding():
    tenth = 0.1 / 60.0 / 60.0  # of an arcsecond

    assert str(Angle(degrees=tenth * -600.75)) == '-0deg 1\' 0.1"'
    assert str(Angle(degrees=tenth * -600.25)) == '-0deg 1\' 0.0"'
    assert str(Angle(degrees=tenth * -599.75)) == '-0deg 1\' 0.0"'
    assert str(Angle(degrees=tenth * -599.25)) == '-0deg 0\' 59.9"'

    assert str(Angle(degrees=tenth * -1.75)) == '-0deg 0\' 0.2"'
    assert str(Angle(degrees=tenth * -1.25)) == '-0deg 0\' 0.1"'
    assert str(Angle(degrees=tenth * -0.75)) == '-0deg 0\' 0.1"'
    assert str(Angle(degrees=tenth * -0.25)) == '-0deg 0\' 0.0"'

    assert str(Angle(degrees=0.0)) == '0deg 0\' 0.0"'

    assert str(Angle(degrees=tenth * 0.25)) == '0deg 0\' 0.0"'
    assert str(Angle(degrees=tenth * 0.75)) == '0deg 0\' 0.1"'
    assert str(Angle(degrees=tenth * 1.25)) == '0deg 0\' 0.1"'
    assert str(Angle(degrees=tenth * 1.75)) == '0deg 0\' 0.2"'

    assert str(Angle(degrees=tenth * 599.25)) == '0deg 0\' 59.9"'
    assert str(Angle(degrees=tenth * 599.75)) == '0deg 1\' 0.0"'
    assert str(Angle(degrees=tenth * 600.25)) == '0deg 1\' 0.0"'
    assert str(Angle(degrees=tenth * 600.75)) == '0deg 1\' 0.1"'
