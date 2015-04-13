import numpy as np
from assay import assert_raises
from pytz import timezone
from skyfield.constants import DAY_S
from skyfield.timelib import JulianDate, takes_julian_date, utc
from datetime import datetime

one_second = 1.0 / DAY_S
epsilon = one_second * 42.0e-6  # 20.1e-6 is theoretical best precision

time_parameter = ['tai', 'tt', 'tdb']
time_value = [(1973, 1, 18, 1, 35, 37.5), 2441700.56640625]

def test_tuple_error_raised():
    def fake_function():
        """ The decorator requires a docstring. """
    decorated_function = takes_julian_date(fake_function)
    tuple_argument = (2014, 1, 23)
    self = None
    with assert_raises(ValueError) as info:
        decorated_function(self, tuple_argument)
    assert 'Are you trying to pass in a tuple' in str(info.exception)

def test_JulianDate_init(time_parameter, time_value):
    kw = {time_parameter: time_value}
    jd = JulianDate(**kw)
    assert getattr(jd, time_parameter) == 2441700.56640625

def test_building_JulianDate_from_utc_tuple_with_array_inside():
    seconds = np.arange(48.0, 58.0, 1.0)
    jd = JulianDate(utc=(1973, 12, 29, 23, 59, seconds))
    assert seconds.shape == jd.shape
    for i, second in enumerate(seconds):
        assert jd.tai[i] == JulianDate(utc=(1973, 12, 29, 23, 59, second)).tai

def test_building_JulianDate_from_naive_datetime_raises_exception():
    with assert_raises(ValueError) as info:
        JulianDate(utc=datetime(1973, 12, 29, 23, 59, 48))
    assert 'import timezone' in str(info.exception)

def test_building_JulianDate_from_single_utc_datetime():
    jd = JulianDate(utc=datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc))
    assert jd.tai == 2442046.5

def test_building_JulianDate_from_list_of_utc_datetimes():
    jd = JulianDate(utc=[
        datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 30, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 31, 23, 59, 48, tzinfo=utc),
        datetime(1974, 1, 1, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 2, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 3, 23, 59, 47, tzinfo=utc),
        ])
    assert (jd.tai == [
        2442046.5, 2442047.5, 2442048.5, 2442049.5, 2442050.5, 2442051.5,
        ]).all()

def test_indexing_julian_date():
    jd = JulianDate(utc=(1974, 10, range(1, 6)))
    assert jd.shape == (5,)
    jd0 = jd[0]
    assert jd.tai[0] == jd0.tai
    assert jd.tt[0] == jd0.tt
    assert jd.tdb[0] == jd0.tdb
    assert jd.ut1[0] == jd0.ut1
    assert jd.delta_t == jd0.delta_t

def test_slicing_julian_date():
    jd = JulianDate(utc=(1974, 10, range(1, 6)))
    assert jd.shape == (5,)
    jd24 = jd[2:4]
    assert jd24.shape == (2,)
    assert (jd.tai[2:4] == jd24.tai).all()
    assert (jd.tt[2:4] == jd24.tt).all()
    assert (jd.tdb[2:4] == jd24.tdb).all()
    assert (jd.ut1[2:4] == jd24.ut1).all()
    assert jd.delta_t == jd24.delta_t

def test_early_utc():
    jd = JulianDate(utc=(1915, 12, 2, 3, 4, 5.6786786))
    assert abs(jd.tt - 2420833.6283317441) < epsilon
    assert jd.utc_iso() == '1915-12-02T03:04:06Z'

def test_astimezone():
    jd = JulianDate(utc=(1969, 7, 20, 20, 18))
    tz = timezone('US/Eastern')
    dt = jd.astimezone(tz)
    assert dt == tz.localize(datetime(1969, 7, 20, 16, 18, 0, 0))

def test_astimezone_and_leap_second():
    jd = JulianDate(utc=(1969, 7, 20, 20, 18))
    tz = timezone('US/Eastern')
    dt, leap_second = jd.astimezone_and_leap_second(tz)
    assert dt == tz.localize(datetime(1969, 7, 20, 16, 18, 0, 0))
    assert leap_second == 0

def test_utc_datetime():
    jd = JulianDate(utc=(1969, 7, 20, 20, 18))
    dt = jd.utc_datetime()
    assert dt == datetime(1969, 7, 20, 20, 18, 0, 0, utc)

def test_utc_datetime_and_leap_second():
    jd = JulianDate(utc=(1969, 7, 20, 20, 18))
    dt, leap_second = jd.utc_datetime_and_leap_second()
    assert dt == datetime(1969, 7, 20, 20, 18, 0, 0, utc)
    assert leap_second == 0

def test_iso_of_decimal_that_rounds_up():
    jd = JulianDate(utc=(1915, 12, 2, 3, 4, 5.6786786))
    assert jd.utc_iso(places=0) == '1915-12-02T03:04:06Z'
    assert jd.utc_iso(places=1) == '1915-12-02T03:04:05.7Z'
    assert jd.utc_iso(places=2) == '1915-12-02T03:04:05.68Z'
    assert jd.utc_iso(places=3) == '1915-12-02T03:04:05.679Z'
    assert jd.utc_iso(places=4) == '1915-12-02T03:04:05.6787Z'

def test_iso_of_decimal_that_rounds_down():
    jd = JulianDate(utc=(2014, 12, 21, 6, 3, 1.234234))
    assert jd.utc_iso(places=0) == '2014-12-21T06:03:01Z'
    assert jd.utc_iso(places=1) == '2014-12-21T06:03:01.2Z'
    assert jd.utc_iso(places=2) == '2014-12-21T06:03:01.23Z'
    assert jd.utc_iso(places=3) == '2014-12-21T06:03:01.234Z'
    assert jd.utc_iso(places=4) == '2014-12-21T06:03:01.2342Z'

def test_iso_of_leap_second_with_fraction():
    jd = JulianDate(utc=(1973, 12, 31, 23, 59, 60.12349))
    assert jd.utc_iso(places=0) == '1973-12-31T23:59:60Z'
    assert jd.utc_iso(places=1) == '1973-12-31T23:59:60.1Z'
    assert jd.utc_iso(places=2) == '1973-12-31T23:59:60.12Z'
    assert jd.utc_iso(places=3) == '1973-12-31T23:59:60.123Z'
    assert jd.utc_iso(places=4) == '1973-12-31T23:59:60.1235Z'

def test_iso_of_array_showing_whole_seconds():
    jd = JulianDate(utc=(1973, 12, 31, 23, 59, np.arange(58.75, 63.1, 0.5)))
    assert jd.utc_iso(places=0) == [
        '1973-12-31T23:59:59Z',
        '1973-12-31T23:59:59Z',
        '1973-12-31T23:59:60Z',
        '1973-12-31T23:59:60Z',
        '1974-01-01T00:00:00Z',
        '1974-01-01T00:00:00Z',
        '1974-01-01T00:00:01Z',
        '1974-01-01T00:00:01Z',
        '1974-01-01T00:00:02Z',
        ]

def test_iso_of_array_showing_fractions():
    jd = JulianDate(utc=(1973, 12, 31, 23, 59, np.arange(58.75, 63.1, 0.5)))
    assert jd.utc_iso(places=2) == [
        '1973-12-31T23:59:58.75Z',
        '1973-12-31T23:59:59.25Z',
        '1973-12-31T23:59:59.75Z',
        '1973-12-31T23:59:60.25Z',
        '1973-12-31T23:59:60.75Z',
        '1974-01-01T00:00:00.25Z',
        '1974-01-01T00:00:00.75Z',
        '1974-01-01T00:00:01.25Z',
        '1974-01-01T00:00:01.75Z',
        ]

def test_jpl_format():
    jd = JulianDate(utc=(range(-300, 301, 100), 7, 1))
    assert jd.utc_jpl() == [
        'B.C. 0301-Jul-01 00:00:00.0000 UT',
        'B.C. 0201-Jul-01 00:00:00.0000 UT',
        'B.C. 0101-Jul-01 00:00:00.0000 UT',
        'B.C. 0001-Jul-01 00:00:00.0000 UT',
        'A.D. 0100-Jul-01 00:00:00.0000 UT',
        'A.D. 0200-Jul-01 00:00:00.0000 UT',
        'A.D. 0300-Jul-01 00:00:00.0000 UT',
        ]

def test_stftime_of_single_date():
    jd = JulianDate(utc=(1973, 12, 31, 23, 59, 60))
    assert jd.utc_strftime('%Y %m %d %H %M %S') == '1973 12 31 23 59 60'

def test_stftime_of_date_array():
    jd = JulianDate(utc=(1973, 12, 31, 23, 59, np.arange(59.0, 61.1, 1.0)))
    assert jd.utc_strftime('%Y %m %d %H %M %S') == [
        '1973 12 31 23 59 59',
        '1973 12 31 23 59 60',
        '1974 01 01 00 00 00',
        ]

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
