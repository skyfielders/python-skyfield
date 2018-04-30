import numpy as np
from skyfield.units import Angle

def test_degree_rounding():
    tenth = 0.1 / 60.0 / 60.0  # of an arcsecond

    assert str(Angle(degrees=tenth * -600.75)) == '-00deg 01\' 00.1"'
    assert str(Angle(degrees=tenth * -600.25)) == '-00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * -599.75)) == '-00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * -599.25)) == '-00deg 00\' 59.9"'

    assert str(Angle(degrees=tenth * -1.75)) == '-00deg 00\' 00.2"'
    assert str(Angle(degrees=tenth * -1.25)) == '-00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * -0.75)) == '-00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * -0.25)) == '-00deg 00\' 00.0"'

    assert str(Angle(degrees=0.0)) == '00deg 00\' 00.0"'

    assert str(Angle(degrees=tenth * 0.25)) == '00deg 00\' 00.0"'
    assert str(Angle(degrees=tenth * 0.75)) == '00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * 1.25)) == '00deg 00\' 00.1"'
    assert str(Angle(degrees=tenth * 1.75)) == '00deg 00\' 00.2"'

    assert str(Angle(degrees=tenth * 599.25)) == '00deg 00\' 59.9"'
    assert str(Angle(degrees=tenth * 599.75)) == '00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * 600.25)) == '00deg 01\' 00.0"'
    assert str(Angle(degrees=tenth * 600.75)) == '00deg 01\' 00.1"'

def test_angle_scalar_strs():
    assert str(Angle(degrees=np.array(91))) == '''91deg 00' 00.0"'''
    assert str(Angle(hours=np.array(12))) == '''12h 00m 00.00s'''

def test_angle_array_strs():
    assert str(Angle(degrees=np.array([90, 91, 92]))) == (
        '''3 values from 90deg 00' 00.0" to 92deg 00' 00.0"'''
        )
    assert str(Angle(hours=np.array([11, 12, 13]))) == (
        '''3 values from 11h 00m 00.00s to 13h 00m 00.00s'''
        )

def test_angle_sexagesimal_args():
    assert str(Angle(degrees=(90,))) == '''90deg 00' 00.0"'''
    assert str(Angle(hours=(12,))) == '''12h 00m 00.00s'''

    assert str(Angle(degrees=(90, 15))) == '''90deg 15' 00.0"'''
    assert str(Angle(hours=(12, 30))) == '''12h 30m 00.00s'''

    assert str(Angle(degrees=(90, 15, 30))) == '''90deg 15' 30.0"'''
    assert str(Angle(hours=(12, 30, 15))) == '''12h 30m 15.00s'''
