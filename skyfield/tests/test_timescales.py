from numpy import searchsorted
from skyfield.constants import DAY_S
from skyfield.io import Cache
from skyfield.timescales import download_leapseconds, julian_date

def test_leap_second():
    cache = Cache('.', days_old=9999)
    leap_dates, leap_offsets = cache.run(download_leapseconds)

    def from_utc_to_tai(year, month, day, hour, minute, second):
        j = julian_date(year, month, day, hour, minute, 0.0)
        i = searchsorted(leap_dates, j, 'right')
        return j + (second + leap_offsets[i]) / DAY_S

    # During 1973 the offset between UTC and TAI was 12.0 seconds, so
    # TAI should reach the first moment of 1974 while the UTC clock is
    # still reading 12s before midnight (60 - 12 = 48).  Happily, the
    # fraction 0.5 can be precisely represented in floating point, so we
    # can use a bare `==` in this assert:

    t0 = from_utc_to_tai(1973, 12, 31, 23, 59, 48.0)
    assert t0 == 2442048.5

    # Here are some more interesting values:

    t1 = from_utc_to_tai(1973, 12, 31, 23, 59, 58.0)
    t2 = from_utc_to_tai(1973, 12, 31, 23, 59, 59.0)
    t3 = from_utc_to_tai(1973, 12, 31, 23, 59, 60.0)
    t4 = from_utc_to_tai(1974, 1, 1, 0, 0, 0.0)
    t5 = from_utc_to_tai(1974, 1, 1, 0, 0, 1.0)

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
