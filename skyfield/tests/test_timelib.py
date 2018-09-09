import numpy as np
from assay import assert_raises
from pytz import timezone
from skyfield import api
from skyfield.constants import DAY_S
from skyfield.timelib import utc, calendar_tuple, julian_date
from datetime import datetime

one_second = 1.0 / DAY_S
epsilon = one_second * 42.0e-6  # 20.1e-6 is theoretical best precision

time_parameter = ['tai', 'tt', 'tdb', 'ut1']
time_value = [(1973, 1, 18, 1, 35, 37.5), 2441700.56640625]

def ts():
    yield api.load.timescale()

def test_time_creation_methods(ts, time_parameter, time_value):
    method = getattr(ts, time_parameter)
    if isinstance(time_value, tuple):
        t = method(*time_value)
    else:
        t = method(jd=time_value) # TODO: deprecate
    assert getattr(t, time_parameter) == 2441700.56640625

def test_time_creation_from_julian_day(ts, time_parameter):
    # Test methods like "tt_jd()" and "tai_jd()".
    method = getattr(ts, time_parameter + '_jd')
    t = method(2441700.56640625)
    assert getattr(t, time_parameter) == 2441700.56640625

time_scale_name = ['utc', 'tai', 'tt', 'tdb']
time_params_with_array = [
    ((2018, 2019, 2020), 3, 25, 13, 1, 10),
    (2018, (3, 4, 5), 25, 13, 1, 10),
    (2018, 3, (25, 26, 27), 13, 1, 10),
    (2018, 3, 25, (13, 14, 15), 1, 10),
    (2018, 3, 25, 13, (1, 2, 3), 10),
    (2018, 3, 25, 13, 1, (10, 11, 12)),
]

def test_time_creation_with_arrays(time_scale_name, time_params_with_array):
    ts = api.load.timescale()
    getattr(ts, time_scale_name)(*time_params_with_array)

def test_timescale_utc_method_with_array_inside(ts):
    seconds = np.arange(48.0, 58.0, 1.0)
    t = ts.utc(1973, 12, 29, 23, 59, seconds)
    assert seconds.shape == t.shape
    for i, second in enumerate(seconds):
        assert t.tai[i] == ts.utc(1973, 12, 29, 23, 59, second).tai

def test_that_building_time_from_naive_datetime_raises_exception(ts):
    with assert_raises(ValueError) as info:
        ts.utc(datetime(1973, 12, 29, 23, 59, 48))
    assert 'import timezone' in str(info.exception)

def test_building_time_from_single_utc_datetime(ts):
    t = ts.utc(datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc))
    assert t.tai == 2442046.5

def test_building_time_from_list_of_utc_datetimes(ts):
    t = ts.utc([
        datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 30, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 31, 23, 59, 48, tzinfo=utc),
        datetime(1974, 1, 1, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 2, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 3, 23, 59, 47, tzinfo=utc),
        ])
    assert (t.tai == [
        2442046.5, 2442047.5, 2442048.5, 2442049.5, 2442050.5, 2442051.5,
        ]).all()

def test_converting_ut1_to_tt(ts):
    ten_thousand_years = 365 * 10000

    jd = api.T0 - ten_thousand_years
    t = ts.ut1(jd=jd)
    del t.ut1                   # force re-computation of UT1
    print(jd - t.ut1)
    assert abs(jd - t.ut1) < 1e-10

    jd = api.T0 + ten_thousand_years
    t = ts.ut1(jd=jd)
    del t.ut1                   # force re-computation of UT1
    print(jd - t.ut1)
    assert abs(jd - t.ut1) < 1e-10

def test_indexing_time(ts):
    t = ts.utc(1974, 10, range(1, 6))
    assert t.shape == (5,)
    t0 = t[0]
    assert t.tai[0] == t0.tai
    assert t.tt[0] == t0.tt
    assert t.tdb[0] == t0.tdb
    assert t.ut1[0] == t0.ut1
    assert t.delta_t[0] == t0.delta_t

def test_slicing_time(ts):
    t = ts.utc(1974, 10, range(1, 6))
    assert t.shape == (5,)
    t24 = t[2:4]
    assert t24.shape == (2,)
    assert (t.tai[2:4] == t24.tai).all()
    assert (t.tt[2:4] == t24.tt).all()
    assert (t.tdb[2:4] == t24.tdb).all()
    assert (t.ut1[2:4] == t24.ut1).all()
    assert (t.delta_t[2:4] == t24.delta_t).all()

def test_early_utc(ts):
    t = ts.utc(1915, 12, 2, 3, 4, 5.6786786)
    assert abs(t.tt - 2420833.6283317441) < epsilon
    assert t.utc_iso() == '1915-12-02T03:04:06Z'

def test_astimezone(ts):
    t = ts.utc(1969, 7, 20, 20, 18)
    tz = timezone('US/Eastern')
    dt = t.astimezone(tz)
    assert dt == tz.localize(datetime(1969, 7, 20, 16, 18, 0, 0))

def test_astimezone_and_leap_second(ts):
    t = ts.utc(1969, 7, 20, 20, 18)
    tz = timezone('US/Eastern')
    dt, leap_second = t.astimezone_and_leap_second(tz)
    assert dt == tz.localize(datetime(1969, 7, 20, 16, 18, 0, 0))
    assert leap_second == 0

def test_utc_datetime(ts):
    t = ts.utc(1969, 7, 20, 20, 18)
    dt = t.utc_datetime()
    assert dt == datetime(1969, 7, 20, 20, 18, 0, 0, utc)

def test_utc_datetime_and_leap_second(ts):
    t = ts.utc(1969, 7, 20, 20, 18)
    dt, leap_second = t.utc_datetime_and_leap_second()
    assert dt == datetime(1969, 7, 20, 20, 18, 0, 0, utc)
    assert leap_second == 0

def test_iso_of_decimal_that_rounds_up(ts):
    t = ts.utc(1915, 12, 2, 3, 4, 5.6786786)
    assert t.utc_iso(places=0) == '1915-12-02T03:04:06Z'
    assert t.utc_iso(places=1) == '1915-12-02T03:04:05.7Z'
    assert t.utc_iso(places=2) == '1915-12-02T03:04:05.68Z'
    assert t.utc_iso(places=3) == '1915-12-02T03:04:05.679Z'
    assert t.utc_iso(places=4) == '1915-12-02T03:04:05.6787Z'

def test_iso_of_decimal_that_rounds_down(ts):
    t = ts.utc(2014, 12, 21, 6, 3, 1.234234)
    assert t.utc_iso(places=0) == '2014-12-21T06:03:01Z'
    assert t.utc_iso(places=1) == '2014-12-21T06:03:01.2Z'
    assert t.utc_iso(places=2) == '2014-12-21T06:03:01.23Z'
    assert t.utc_iso(places=3) == '2014-12-21T06:03:01.234Z'
    assert t.utc_iso(places=4) == '2014-12-21T06:03:01.2342Z'

def test_iso_of_leap_second_with_fraction(ts):
    t = ts.utc(1973, 12, 31, 23, 59, 60.12349)
    assert t.utc_iso(places=0) == '1973-12-31T23:59:60Z'
    assert t.utc_iso(places=1) == '1973-12-31T23:59:60.1Z'
    assert t.utc_iso(places=2) == '1973-12-31T23:59:60.12Z'
    assert t.utc_iso(places=3) == '1973-12-31T23:59:60.123Z'
    assert t.utc_iso(places=4) == '1973-12-31T23:59:60.1235Z'

def test_iso_of_array_showing_whole_seconds(ts):
    t = ts.utc(1973, 12, 31, 23, 59, np.arange(58.75, 63.1, 0.5))
    assert t.utc_iso(places=0) == [
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

def test_iso_of_array_showing_fractions(ts):
    t = ts.utc(1973, 12, 31, 23, 59, np.arange(58.75, 63.1, 0.5))
    assert t.utc_iso(places=2) == [
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

def test_jpl_format(ts):
    t = ts.utc(range(-300, 301, 100), 7, 1)
    assert t.utc_jpl() == [
        'B.C. 0301-Jul-01 00:00:00.0000 UT',
        'B.C. 0201-Jul-01 00:00:00.0000 UT',
        'B.C. 0101-Jul-01 00:00:00.0000 UT',
        'B.C. 0001-Jul-01 00:00:00.0000 UT',
        'A.D. 0100-Jul-01 00:00:00.0000 UT',
        'A.D. 0200-Jul-01 00:00:00.0000 UT',
        'A.D. 0300-Jul-01 00:00:00.0000 UT',
        ]

def test_stftime_of_single_date(ts):
    t = ts.utc(1973, 12, 31, 23, 59, 60)
    assert t.utc_strftime('%Y %m %d %H %M %S') == '1973 12 31 23 59 60'

def test_stftime_of_date_array(ts):
    t = ts.utc(1973, 12, 31, 23, 59, np.arange(59.0, 61.1, 1.0))
    assert t.utc_strftime('%Y %m %d %H %M %S') == [
        '1973 12 31 23 59 59',
        '1973 12 31 23 59 60',
        '1974 01 01 00 00 00',
        ]

def test_leap_second(ts):

    # During 1973 the offset between UTC and TAI was 12.0 seconds, so
    # TAI should reach the first moment of 1974 while the UTC clock is
    # still reading 12s before midnight (60 - 12 = 48).  Happily, the
    # fraction 0.5 can be precisely represented in floating point, so we
    # can use a bare `==` in this assert:

    t0 = ts.utc(1973, 12, 31, 23, 59, 48.0).tai
    assert t0 == 2442048.5

    # Here are some more interesting values:

    t1 = ts.utc(1973, 12, 31, 23, 59, 58.0).tai
    t2 = ts.utc(1973, 12, 31, 23, 59, 59.0).tai
    t3 = ts.utc(1973, 12, 31, 23, 59, 60.0).tai
    t4 = ts.utc(1974, 1, 1, 0, 0, 0.0).tai
    t5 = ts.utc(1974, 1, 1, 0, 0, 1.0).tai

    # The step from 23:59:59 to 0:00:00 is here a two-second step,
    # because of the leap second 23:59:60 that falls in between:

    assert abs(t4 - t2 - 2.0 * one_second) < epsilon

    # Thus, the five dates given above are all one second apart:

    assert abs(t2 - t1 - one_second) < epsilon
    assert abs(t3 - t2 - one_second) < epsilon
    assert abs(t4 - t3 - one_second) < epsilon
    assert abs(t5 - t4 - one_second) < epsilon

    # And all these dates can be converted back to UTC.

    assert ts.tai(jd=t0).utc_iso() == '1973-12-31T23:59:48Z'
    assert ts.tai(jd=t1).utc_iso() == '1973-12-31T23:59:58Z'
    assert ts.tai(jd=t2).utc_iso() == '1973-12-31T23:59:59Z'
    assert ts.tai(jd=t3).utc_iso() == '1973-12-31T23:59:60Z'
    assert ts.tai(jd=t4).utc_iso() == '1974-01-01T00:00:00Z'
    assert ts.tai(jd=t5).utc_iso() == '1974-01-01T00:00:01Z'

def test_delta_t(ts):
    # Check delta_t calculation around year 2000/1/1 (from IERS tables this is 63.8285)
    t = ts.utc(2000, 1, 1, 0, 0, 0)
    assert abs(t.delta_t - 63.8285) < 1e-5

    # Check historic value. Compare to the table in Morrison and
    # Stephenson 2004, the tolerance is 2 sigma
    t = ts.utc(year=1000)
    assert abs(t.delta_t - 1570.0) < 110.0

    # Check future value. Should be calculated by Morrison and
    # Stephenson formula. For 2320 (t=5 cy) should be: -20 + 32 * 5**2
    t = ts.utc(year=2320)
    assert abs(t.delta_t + 20.0 - (32.0 * 5.0**2)) < 1.0

def test_J(ts):
    assert ts.tt(2000, 1, 1.5).J == 2000.0
    assert ts.tt(1900, 1, 0.5).J == 1900.0

def test_time_repr(ts):

    # Check that repr return is a str (this is required on Python 2,
    # unicode is not allowed)
    assert isinstance(repr(ts.utc(year=2000)), str)

    # Check array conversion
    assert isinstance(repr(ts.utc(year=range(2000, 2010))), str)

    assert repr(ts.tt_jd(1)) == '<Time tt=1.0>'
    assert repr(ts.tt_jd([])) == '<Time tt=[]>'
    assert repr(ts.tt_jd([1])) == '<Time tt=[1]>'
    assert repr(ts.tt_jd([1, 2])) == '<Time tt=[1 2]>'
    assert repr(ts.tt_jd([1, 2, 3])) == '<Time tt=[1 2 3]>'
    assert repr(ts.tt_jd([1, 2, 3, 4])) == '<Time tt=[1 ... 4] len=4>'

def test_jd_calendar():

    import numbers

    # Check a specific instance (using UNIX epoch here, though that's
    # an arbitrary choice)
    jd_unix = 2440587.5
    cal_unix = (1970, 1, 1, 0, 0, 0.0)

    cal = calendar_tuple(jd_unix)

    # Check that the returned value is correct
    assert cal == cal_unix

    # Check that all the return types are correct
    assert isinstance(cal[0], numbers.Integral)  # Year
    assert isinstance(cal[1], numbers.Integral)  # Month
    assert isinstance(cal[2], numbers.Integral)  # Day
    assert isinstance(cal[3], numbers.Integral)  # Hour
    assert isinstance(cal[4], numbers.Integral)  # Minute
    assert isinstance(cal[5], numbers.Real)  # Second

    # Check backward conversion
    assert julian_date(*cal) == jd_unix

    # Check array conversion components
    jd_array = jd_unix + np.arange(5.0)
    cal_array = calendar_tuple(jd_array)
    
    assert (cal_array[0] == 1970).all()
    assert (cal_array[1] == 1).all()
    assert (cal_array[2] == np.arange(1, 6)).all()
    assert (cal_array[3] == 0).all()
    assert (cal_array[4] == 0).all()
    assert (cal_array[5] == 0.0).all()

    # Check reversal of array
    assert (julian_date(*cal_array) == jd_array).all()
