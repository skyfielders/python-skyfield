from skyfield.constants import DAY_S
from skyfield.timescales import JulianDate

one_second = 1.0 / DAY_S
epsilon = one_second * 42.0e-6  # 20.1e-6 is theoretical best precision

def test_early_utc():
    jd = JulianDate(utc=(1915, 12, 2, 3, 4, 5.6786786))
    assert abs(jd.tt - 2420833.6283317441) < epsilon
    assert jd.utc_iso() == '1915-12-02T03:04:06Z'

def test_iso_decimal_that_rounds_up():
    jd = JulianDate(utc=(1915, 12, 2, 3, 4, 5.6786786))
    assert jd.utc_iso(places=0) == '1915-12-02T03:04:06Z'
    assert jd.utc_iso(places=1) == '1915-12-02T03:04:05.7Z'
    assert jd.utc_iso(places=2) == '1915-12-02T03:04:05.68Z'
    assert jd.utc_iso(places=3) == '1915-12-02T03:04:05.679Z'
    assert jd.utc_iso(places=4) == '1915-12-02T03:04:05.6787Z'

def test_iso_decimal_that_rounds_down():
    jd = JulianDate(utc=(2014, 12, 21, 6, 3, 1.234234))
    assert jd.utc_iso(places=0) == '2014-12-21T06:03:01Z'
    assert jd.utc_iso(places=1) == '2014-12-21T06:03:01.2Z'
    assert jd.utc_iso(places=2) == '2014-12-21T06:03:01.23Z'
    assert jd.utc_iso(places=3) == '2014-12-21T06:03:01.234Z'
    assert jd.utc_iso(places=4) == '2014-12-21T06:03:01.2342Z'

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

    assert abs(t4 - t2 - 2.0 * one_second) < epsilon

    # Thus, the five dates given above are all one second apart:

    assert abs(t2 - t1 - one_second) < epsilon
    assert abs(t3 - t2 - one_second) < epsilon
    assert abs(t4 - t3 - one_second) < epsilon
    assert abs(t5 - t4 - one_second) < epsilon

    # And all these dates can be converted back to UTC.

    assert JulianDate(tai=t0).utc_iso() == '1973-12-31T23:59:48Z'
    assert JulianDate(tai=t1).utc_iso() == '1973-12-31T23:59:58Z'
    assert JulianDate(tai=t2).utc_iso() == '1973-12-31T23:59:59Z'
    assert JulianDate(tai=t3).utc_iso() == '1973-12-31T23:59:60Z'
    assert JulianDate(tai=t4).utc_iso() == '1974-01-01T00:00:00Z'
    assert JulianDate(tai=t5).utc_iso() == '1974-01-01T00:00:01Z'
