import datetime as dt_module
import numbers
import numpy as np
import sys
from assay import assert_raises
from pytz import timezone
from skyfield import api
from skyfield.constants import DAY_S, T0
from skyfield.timelib import (
    GREGORIAN_START, GREGORIAN_START_ENGLAND, Time,
    calendar_tuple, compute_calendar_date, julian_date, julian_day, utc,
)
from datetime import datetime

one_second = 1.0 / DAY_S
epsilon = one_second * 42.0e-6  # 20.1e-6 is theoretical best precision

continuous_timescale = ['tai', 'tt', 'tdb', 'ut1']
time_scale_name = ['utc', 'tai', 'tt', 'tdb', 'ut1']
time_value = [(1973, 1, 18, 1, 35, 37.5), 2441700.56640625]

def ts():
    yield api.load.timescale()

def a(*args):
    return np.array(args)

def all_kinds_of_time_array(ts):
    yield ts.utc(2020, 10, [8, 9])
    yield ts.tai(2020, 10, [8, 9])
    yield ts.tt(2020, 10, [8, 9])
    yield ts.tdb(2020, 10, [8, 9])
    yield ts.ut1(2020, 10, [8, 9])

    jd = a(2459130.5, 2459131.5)

    yield ts.tai_jd(jd)
    yield ts.tt_jd(jd)
    yield ts.tdb_jd(jd)
    yield ts.ut1_jd(jd)
    yield Time(ts, jd)

    for jd, fraction in [
            (a(2459130.5, 2459131.5), 0.25),
            (2459130.5, a(0.0, 0.25)),
            (a(2459130.5, 2459131.5), a(0.0, 0.25)),
    ]:
        yield ts.tai_jd(jd, fraction)
        yield ts.tt_jd(jd, fraction)
        yield ts.tdb_jd(jd, fraction)
        # yield ts.ut1_jd(jd, fraction)  # not yet supported

    # We only support direct Time instantiation for the final case,
    # where jd and fraction are arrays that already agree in their
    # dimensions.

    yield Time(ts, jd, fraction)

def test_time_creation_methods(ts, continuous_timescale, time_value):
    method = getattr(ts, continuous_timescale)
    if isinstance(time_value, tuple):
        t = method(*time_value)
    else:
        t = method(jd=time_value) # TODO: deprecate
    assert getattr(t, continuous_timescale) == 2441700.56640625

    # Also go ahead and test the calendar and formatting operations.

    tup = getattr(t, continuous_timescale + '_calendar')()
    assert tup == (1973, 1, 18, 1, 35, 37.5)

    strftime = getattr(t, continuous_timescale + '_strftime')
    string = strftime()
    assert string == '1973-01-18 01:35:38 ' + continuous_timescale.upper()

    if sys.version_info <= (3,):
        return  # we do not currently support %f under Python 2

    string = strftime('%S.%f')
    assert string == '37.500000'

def test_months_overflow_correctly(ts):
    assert ts.tt(2020, -1).tt_strftime('%Y-%m') == '2019-11'
    assert ts.tt(2020, 15).tt_strftime('%Y-%m') == '2021-03'
    assert ts.tt(2020, [-1, 0, 1, 13, 14, 15]).tt_strftime('%Y-%m') == [
        '2019-11', '2019-12', '2020-01', '2021-01', '2021-02', '2021-03',
    ]

def test_days_overflow_correctly(ts):
    months = range(1, 13)
    assert ts.tt(2020, months, -1).tt_strftime('%Y-%m-%d') == [
        '2019-12-30', '2020-01-30', '2020-02-28', '2020-03-30',
        '2020-04-29', '2020-05-30', '2020-06-29', '2020-07-30',
        '2020-08-30', '2020-09-29', '2020-10-30', '2020-11-29',
    ]
    assert ts.tt(2020, months, 0).tt_strftime('%Y-%m-%d') == [
        '2019-12-31', '2020-01-31', '2020-02-29', '2020-03-31',
        '2020-04-30', '2020-05-31', '2020-06-30', '2020-07-31',
        '2020-08-31', '2020-09-30', '2020-10-31', '2020-11-30',
    ]
    assert ts.tt(2020, months, 32).tt_strftime('%Y-%m-%d') == [
        '2020-02-01', '2020-03-03', '2020-04-01', '2020-05-02',
        '2020-06-01', '2020-07-02', '2020-08-01', '2020-09-01',
        '2020-10-02', '2020-11-01', '2020-12-02', '2021-01-01',
    ]

def test_time_can_be_indexed(ts):
    for t in all_kinds_of_time_array(ts):
        t[0]

def test_is_time_iterable(ts, time_scale_name):
    t = getattr(ts, time_scale_name)(2020, 9, (25, 26))
    for item in t:
        pass

def test_strftime_on_prehistoric_dates(ts):
    if sys.version_info <= (3,):
        return  # Python 2 time.strftime() complains about negative years

    t = ts.tt(-746, 2, 26)
    assert t.utc_strftime('%Y %S') == '-746 18'
    assert t.ut1_strftime('%Y %S') == '-746 14'
    assert t.tai_strftime('%Y %S') == '-746 28'
    assert t.tt_strftime('%Y %S') == '-746 00'
    assert t.tdb_strftime('%Y %S') == '-746 00'

    t = ts.tt(-746, 2, [26, 26])
    assert t.utc_strftime('%Y %S') == ['-746 18'] * 2
    assert t.ut1_strftime('%Y %S') == ['-746 14'] * 2
    assert t.tai_strftime('%Y %S') == ['-746 28'] * 2
    assert t.tt_strftime('%Y %S') == ['-746 00'] * 2
    assert t.tdb_strftime('%Y %S') == ['-746 00'] * 2

def test_strftime_with_microseconds(ts):
    if sys.version_info <= (3,):
        return  # we do not currently support %f under Python 2

    t = ts.tt(2020, 9, 12)
    assert t.utc_strftime('%Y %S %f') == '2020 50 816000'
    assert t.ut1_strftime('%Y %S %f') == '2020 49 776521'
    assert t.tai_strftime('%Y %S %f') == '2020 27 816000'
    assert t.tt_strftime('%Y %S %f') == '2020 00 000000'
    assert t.tdb_strftime('%Y %S %f') == '2020 59 998446'

    t = ts.tt(2020, 9, [12, 12])
    assert t.utc_strftime('%Y %S %f') == ['2020 50 816000'] * 2
    assert t.ut1_strftime('%Y %S %f') == ['2020 49 776521'] * 2
    assert t.tai_strftime('%Y %S %f') == ['2020 27 816000'] * 2
    assert t.tt_strftime('%Y %S %f') == ['2020 00 000000'] * 2
    assert t.tdb_strftime('%Y %S %f') == ['2020 59 998446'] * 2

def test_tai_fraction_loses_no_precision(ts):
    t = ts.tai_jd(2459008.0, 0.0123456789)
    assert t.whole == 2459008.0
    assert t.tai_fraction == 0.0123456789

def test_tdb_fraction_loses_no_precision(ts):
    t = ts.tdb_jd(2459008.0, 0.0123456789)
    assert t.whole == 2459008.0
    assert t.tdb_fraction == 0.0123456789

def test_tai_seconds_preserve_10_decimal_places_in_calendar_seconds(ts):
    t = ts.tai(2020, 6, 7, 2, 2, 12.0123456789)
    c = t.tai_calendar()
    assert c[:5] == (2020, 6, 7, 2, 2)
    assert '%.10f' % c[5] == '12.0123456789'

def test_tt_seconds_preserve_10_decimal_places_in_calendar_seconds(ts):
    t = ts.tt(2020, 6, 7, 2, 2, 12.0123456789)
    c = t.tt_calendar()
    assert c[:5] == (2020, 6, 7, 2, 2)
    assert '%.10f' % c[5] == '12.0123456789'

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
    print(time_scale_name)
    t = getattr(ts, time_scale_name)(*time_params_with_array)
    t.utc_jpl()  # a reasonably complicated operation

def test_timescale_utc_method_with_array_inside(ts):
    seconds = np.arange(48.0, 58.0, 1.0)
    t = ts.utc(1973, 12, 29, 23, 59, seconds)
    assert seconds.shape == t.shape
    for i, second in enumerate(seconds):
        assert t.tai[i] == ts.utc(1973, 12, 29, 23, 59, second).tai

def test_that_building_time_from_naive_datetime_raises_exception(ts):
    with assert_raises(ValueError) as info:
        ts.from_datetime(datetime(1973, 12, 29, 23, 59, 48))
    assert 'import timezone' in str(info.exception)

def test_building_time_from_single_utc_datetime(ts):
    t = ts.from_datetime(datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc))
    assert t.tai == 2442046.5
    t = ts.utc(datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc))
    assert t.tai == 2442046.5

def test_building_time_from_single_utc_datetime_with_timezone(ts):
    tz = timezone('US/Eastern')
    t = ts.from_datetime(tz.localize(datetime(2020, 5, 10, 12, 44, 13, 797865)))
    dt, leap_second = t.utc_datetime_and_leap_second()
    assert dt == datetime(2020, 5, 10, 16, 44, 13, 797865, tzinfo=utc)
    assert leap_second == 0

def test_building_time_from_list_of_utc_datetimes(ts):
    datetimes = [
        datetime(1973, 12, 29, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 30, 23, 59, 48, tzinfo=utc),
        datetime(1973, 12, 31, 23, 59, 48, tzinfo=utc),
        datetime(1974, 1, 1, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 2, 23, 59, 47, tzinfo=utc),
        datetime(1974, 1, 3, 23, 59, 47, tzinfo=utc),
    ]
    t = ts.from_datetimes(datetimes)
    assert list(t.tai) == [
        2442046.5, 2442047.5, 2442048.5, 2442049.5, 2442050.5, 2442051.5,
    ]
    t = ts.utc(datetimes)
    assert list(t.tai) == [
        2442046.5, 2442047.5, 2442048.5, 2442049.5, 2442050.5, 2442051.5,
    ]

def test_building_time_from_python_date(ts):
    d = dt_module.date(2020, 7, 22)
    t = ts.utc(d)
    assert t.utc == (2020, 7, 22, 0, 0, 0.0)

def test_converting_ut1_to_tt(ts):
    ten_thousand_years = 365 * 10000

    jd = api.T0 - ten_thousand_years
    t = ts.ut1(jd=jd)
    del t.ut1_fraction          # force re-computation of UT1
    print(jd - t.ut1)
    assert abs(jd - t.ut1) < 1e-10

    jd = api.T0 + ten_thousand_years
    t = ts.ut1(jd=jd)
    del t.ut1_fraction          # force re-computation of UT1
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

def test_toordinal(ts):
    t = ts.utc(1973, 12, 31, 11, 59, 60)
    assert t.toordinal() == 720623.5

def test_utc_datetime(ts):
    t = ts.utc(1969, 7, 20, 20, 18, 42.186479)
    dt = t.utc_datetime()
    assert dt == datetime(1969, 7, 20, 20, 18, 42, 186479, utc)

def test_utc_datetime_and_leap_second(ts):
    t = ts.utc(1969, 7, 20, 20, 18)
    dt, leap_second = t.utc_datetime_and_leap_second()
    assert dt == datetime(1969, 7, 20, 20, 18, 0, 0, utc)
    assert leap_second == 0

def test_utc_datetime_microseconds_round_trip(ts):
    dt = datetime(2020, 5, 10, 11, 50, 9, 727799, tzinfo=utc)
    t = ts.from_datetime(dt)
    dt2, leap_second = t.utc_datetime_and_leap_second()
    assert dt2 == dt
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

def test_strftime_of_a_leap_second(ts):
    t = ts.utc(1973, 12, 31, 23, 59, 60)
    assert t.utc_strftime('%Y %m %d %H %M %S') == '1973 12 31 23 59 60'

def test_strftime_of_date_array_over_a_leap_second(ts):
    t = ts.utc(1973, 12, 31, 23, 59, np.arange(59.0, 61.1, 1.0))
    assert t.utc_strftime('%a %Y %m %d %H %M %S') == [
        'Mon 1973 12 31 23 59 59',
        'Mon 1973 12 31 23 59 60',
        'Tue 1974 01 01 00 00 00',
    ]

def test_strftime_day_of_year(ts):
    # Based on example date at https://strftime.org/
    assert ts.utc(2013, 9, 29).utc_strftime('%j') == '272'
    assert ts.utc(2013, 9, 30).utc_strftime('%j') == '273'
    assert ts.utc(2013, 9, 30, 23, 59).utc_strftime('%j') == '273'
    assert ts.utc(2013, 9, 30, 23, 60).utc_strftime('%j') == '274'

    assert ts.utc(2013, 9, [29, 30]).utc_strftime('%j') == ['272', '273']

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

    # Otherwise, the five dates given above are all one second apart:

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

def test_dut1(ts):
    # Roughly agreeing with tables on NIST website
    t = ts.utc(2017, 3, 30)
    assert str(t.dut1)[:3] == '0.4'

    t = ts.utc(2018, 9, 21)
    assert str(t.dut1)[:3] == '0.0'

    t = ts.utc(2019, 5, 2)
    assert str(t.dut1)[:4] == '-0.1'  # table says -0.2

def test_J(ts):
    assert ts.J(2000).tt == T0
    assert ts.J(1900).tt == T0 - 36525.0
    assert (ts.J([1900, 2000]).tt == [T0 - 36525.0, T0]).all()
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
    assert repr(ts.tt_jd([1])) == '<Time tt=[1.]>'
    assert repr(ts.tt_jd([1, 2])) == '<Time tt=[1. 2.]>'
    assert repr(ts.tt_jd([1, 2, 3])) == '<Time tt=[1. 2. 3.]>'
    assert repr(ts.tt_jd([1, 2, 3, 4])) == '<Time tt=[1.0 ... 4.0] len=4>'

def test_jd_calendar():
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

def test_raw_julian_gregorian_cutover():
    gregory = 2299161
    assert compute_calendar_date(gregory - 2, gregory) == (1582, 10, 3)
    assert compute_calendar_date(gregory - 1, gregory) == (1582, 10, 4)
    assert compute_calendar_date(gregory + 0, gregory) == (1582, 10, 15)
    assert compute_calendar_date(gregory + 1, gregory) == (1582, 10, 16)

    jd = np.arange(gregory - 2, gregory + 2)
    assert [list(a) for a in compute_calendar_date(jd, gregory)] == [
        [1582, 1582, 1582, 1582],
        [10, 10, 10, 10],
        [3, 4, 15, 16],
    ]

    assert julian_day(1582, 10, 3, gregory) == (gregory - 2)
    assert julian_day(1582, 10, 4, gregory) == (gregory - 1)
    assert julian_day(1582, 10, 15, gregory) == (gregory + 0)
    assert julian_day(1582, 10, 16, gregory) == (gregory + 1)

    days = [3, 4, 15, 16]
    assert list(julian_day(1582, 10, np.array(days), gregory)) == [
        2299159, 2299160, 2299161, 2299162,
    ]

def test_constructor_julian_gregorian_cutover(ts, time_scale_name):
    if sys.version_info <= (3,):
        return  # Python 2 time.strftime() complains about the year 1582

    def jd(y, m, d):
        t = getattr(ts, time_scale_name)(y, m, d)
        if time_scale_name == 'utc':
            return sum(t._utc_float(0.0))
        return getattr(t, time_scale_name)

    ts = api.load.timescale()

    assert jd(1582, 10, 4) == 2299149.5
    assert jd(1582, 10, 15) == 2299160.5
    assert jd(1752, 9, 2) == 2361209.5
    assert jd(1752, 9, 14) == 2361221.5

    ts.julian_calendar_cutoff = GREGORIAN_START

    assert jd(1582, 10, 4) == 2299159.5
    assert jd(1582, 10, 15) == 2299160.5
    assert jd(1752, 9, 2) == 2361209.5
    assert jd(1752, 9, 14) == 2361221.5

    ts.julian_calendar_cutoff = GREGORIAN_START_ENGLAND

    assert jd(1582, 10, 4) == 2299159.5
    assert jd(1582, 10, 15) == 2299170.5
    assert jd(1752, 9, 2) == 2361220.5
    assert jd(1752, 9, 14) == 2361221.5

def test_calendar_tuple_julian_gregorian_cutover(ts, time_scale_name):
    if sys.version_info <= (3,):
        return  # Python 2 time.strftime() complains about the year 1582

    def ymd(jd):
        t = ts.tt_jd(jd, 0.1)
        if time_scale_name == 'utc':
            return t.utc[:3]
        return getattr(t, time_scale_name + '_calendar')()[:3]

    ts = api.load.timescale()

    assert ymd(2299149.5) == (1582, 10, 4)
    assert ymd(2299160.5) == (1582, 10, 15)
    assert ymd(2361209.5) == (1752, 9, 2)
    assert ymd(2361221.5) == (1752, 9, 14)

    ts.julian_calendar_cutoff = GREGORIAN_START

    assert ymd(2299159.5) == (1582, 10, 4)
    assert ymd(2299160.5) == (1582, 10, 15)
    assert ymd(2361209.5) == (1752, 9, 2)
    assert ymd(2361221.5) == (1752, 9, 14)

    ts.julian_calendar_cutoff = GREGORIAN_START_ENGLAND

    assert ymd(2299159.5) == (1582, 10, 4)
    assert ymd(2299170.5) == (1582, 10, 15)
    assert ymd(2361220.5) == (1752, 9, 2)
    assert ymd(2361221.5) == (1752, 9, 14)

def test_strftime_julian_gregorian_cutover(ts, time_scale_name):
    if sys.version_info <= (3,):
        return  # Python 2 time.strftime() complains about the year 1582

    def ymd(jd):
        t = ts.tt_jd(jd, 0.1)
        return getattr(t, time_scale_name + '_strftime')('%Y %m %d')

    ts = api.load.timescale()

    assert ymd(2299149.5) == '1582 10 04'
    assert ymd(2299160.5) == '1582 10 15'
    assert ymd(2361209.5) == '1752 09 02'
    assert ymd(2361221.5) == '1752 09 14'

    ts.julian_calendar_cutoff = GREGORIAN_START

    assert ymd(2299159.5) == '1582 10 04'
    assert ymd(2299160.5) == '1582 10 15'
    assert ymd(2361209.5) == '1752 09 02'
    assert ymd(2361221.5) == '1752 09 14'

    ts.julian_calendar_cutoff = GREGORIAN_START_ENGLAND

    assert ymd(2299159.5) == '1582 10 04'
    assert ymd(2299170.5) == '1582 10 15'
    assert ymd(2361220.5) == '1752 09 02'
    assert ymd(2361221.5) == '1752 09 14'

def test_time_equality(ts):
    t0 = ts.tt_jd(2459008.5, 0.125)
    t1 = ts.tt_jd(2459008.0, 0.625)
    assert t0 == t1
    assert t1 - t0 == 0.0
    assert hash(t0) == hash(t1)

    t2 = ts.tt_jd(2459008.0, 0.6251)
    assert t2 != t0
    assert t2 - t0 > 0
    assert hash(t0) != hash(t2)
