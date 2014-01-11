from skyfield.constants import DAY_S
from skyfield.timescales import JulianDate

def test_leap_second():

    # During 1973 the offset between UTC and TAI was 12.0 seconds, so
    # TAI should reach the first moment of 1974 while the UTC clock is
    # still reading 12s before midnight (60 - 12 = 48).  Happily, the
    # fraction 0.5 can be precisely represented in floating point, so we
    # can use a bare `==` in this assert:

    t0 = JulianDate(utc=(1973, 12, 31, 23, 59, 48.0)).tai
    assert t0 == 2442048.5

    # Here are some more interesting values:

    t1 = JulianDate(utc=(1973, 12, 31, 23, 59, 58.0)).tai
    t2 = JulianDate(utc=(1973, 12, 31, 23, 59, 59.0)).tai
    t3 = JulianDate(utc=(1973, 12, 31, 23, 59, 60.0)).tai
    t4 = JulianDate(utc=(1974, 1, 1, 0, 0, 0.0)).tai
    t5 = JulianDate(utc=(1974, 1, 1, 0, 0, 1.0)).tai

    # The step from 23:59:59 to 0:00:00 is here a two-second step,
    # because of the leap second 23:59:60 that falls in between:

    one_second = 1.0 / DAY_S
    epsilon = one_second * 42.0e-6  # 20.1e-6 is theoretical best precision

    assert abs(t4 - t2 - 2.0 * one_second) < epsilon

    # Thus, the five dates given above are all one second apart:

    assert abs(t2 - t1 - one_second) < epsilon
    assert abs(t3 - t2 - one_second) < epsilon
    assert abs(t4 - t3 - one_second) < epsilon
    assert abs(t5 - t4 - one_second) < epsilon

    # And all five Julian dates can be converted back to UTC.

    assert JulianDate(tai=t1).utc_iso() == '1973-12-31T23:59:58Z'
    assert JulianDate(tai=t2).utc_iso() == '1973-12-31T23:59:59Z'
    assert JulianDate(tai=t3).utc_iso() == '1973-12-31T23:59:60Z'
    assert JulianDate(tai=t4).utc_iso() == '1974-01-01T00:00:00Z'
    assert JulianDate(tai=t5).utc_iso() == '1974-01-01T00:00:01Z'
