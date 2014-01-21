from datetime import date, datetime, timedelta, tzinfo
from numpy import array, einsum, rollaxis, searchsorted, sin, zeros_like
from time import strftime
from .constants import T0, DAY_S
from .framelib import ICRS_to_J2000 as B
from .nutationlib import compute_nutation
from .precessionlib import compute_precession

try:
    from pytz import utc
except ImportError:

    class UTC(tzinfo):
        'UTC'
        zero = timedelta(0)
        def utcoffset(self, dt):
            return self.zero
        def tzname(self, dt):
            return 'UTC'
        def dst(self, dt):
            return self.zero

    utc = UTC()

# Much of the following code is adapted from the USNO's "novas.c".

_half_second = 0.5 / DAY_S
_half_millisecond = 0.5e-3 / DAY_S
_half_microsecond = 0.5e-6 / DAY_S
_months = array(['Month zero', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

extra_documentation = """

        This routine takes a date as its argument.  You can either
        provide a `jd=` keyword argument with a `JulianDate` you have
        built yourself, or use one of these keyword arguments::

            # Coordinated Universal Time
            utc=(1973, 12, 29, 23, 59, 48.0)
            utc=datetime(1973, 12, 29, 23, 59, 48.0)

            # International Atomic Time
            tai=2442046.5

            # Terrestrial Time
            tt=2442046.5

"""

def takes_julian_date(function):
    """Wrap `function` so it accepts the standard Julian date arguments.

    A function that takes two arguments, `self` and `jd`, may be wrapped
    with this decorator if it wants to support optional auto-creation of
    its `jd` argument by accepting all of the same keyword arguments
    that the JulianDate constructor itself supports.

    """
    def wrapper(self, jd=None, utc=None, tai=None, tt=None,
                delta_t=0.0, cache=None):
        if jd is None:
            jd = JulianDate(utc=utc, tai=tai, tt=tt,
                            delta_t=delta_t, cache=cache)
        else:
            pass  # TODO: verify that they provided a JulianDate instance
        return function(self, jd)
    wrapper.__name__ = function.__name__
    synopsis, blank_line, description = function.__doc__.partition('\n\n')
    wrapper.__doc__ = ''.join(synopsis + extra_documentation + description)
    return wrapper

def _wrap(a):
    if hasattr(a, 'shape'):
        return a
    if hasattr(a, '__len__'):
        return array(a)
    return array((a,))

def _to_array(a):
    if not hasattr(a, 'shape'):
        a = array(a)
    return a

tt_minus_tai = array(32.184 / DAY_S)

class JulianDate(object):
    """Julian date.

    :param utc: Coordinated Universal Time
    :param tai: International Atomic Time
    :param tt: Terrestrial Time
    :param delta_t: Difference between Terrestrial Time and UT1
    :param cache: Cache for automatically fetching `delta_t` table

    """
    def __init__(self, utc=None, tai=None, tt=None, delta_t=0.0, cache=None):

        self.delta_t = _to_array(delta_t)

        if cache is None:
            from skyfield.data import cache
        self.cache = cache

        if tai is None and utc is not None:
            leap_dates, leap_offsets = cache.run(usno_leapseconds)
            if isinstance(utc, datetime):
                tai = _utc_datetime_to_tai(leap_dates, leap_offsets, utc)
            elif isinstance(utc, date):
                tai = _utc_date_to_tai(leap_dates, leap_offsets, utc)
            elif isinstance(utc, tuple):
                values = [_to_array(value) for value in utc]
                tai = _utc_to_tai(leap_dates, leap_offsets, *values)
            else:
                tai = array([
                    _utc_datetime_to_tai(leap_dates, leap_offsets, dt)
                    for dt in utc])

        if tai is not None:
            self.tai = _to_array(tai)
            if tt is None:
                tt = tai + tt_minus_tai

        if tt is None:
            raise ValueError('you must supply either utc, tai, or tt when'
                             ' building a JulianDate')
        elif isinstance(tt, tuple):
            tt = julian_date(*tt)

        self.tt = _to_array(tt)
        self.shape = self.tt.shape
        self.delta_t = delta_t

    def astimezone(self, tz):
        """Return a ``datetime`` for this value in a ``pytz`` timezone.

        The third-party ``pytz`` package must be installed separately.
        If this value is a date array, then a sequence of datetimes is
        returned.

        """
        dt, leap_second = self.utc_datetime()
        normalize = getattr(tz, 'normalize', lambda d: d)
        if self.shape:
            dt = [normalize(d.astimezone(tz)) for d in dt]
        else:
            dt = normalize(dt.astimezone(tz))
        return dt, leap_second

    def utc_datetime(self):
        """Return a UTC ``datetime`` for this value.

        If the third-party ``pytz`` package is available, then its
        ``utc`` timezone is used.  Otherwise, a Skyfield built-in
        alternative is used; the result should be the same either way.
        If this value is a date array, then a sequence of datetimes is
        returned.

        """
        year, month, day, hour, minute, second = self._utc(_half_millisecond)
        second, fraction = divmod(second, 1.0)
        second = second.astype(int)
        leap_second = second // 60
        second -= leap_second
        milli = (fraction * 1000).astype(int) * 1000
        if self.shape:
            utcs = [utc] * self.shape[0]
            argsets = zip(year, month, day, hour, minute, second, milli, utcs)
            dt = [datetime(*args) for args in argsets]
        else:
            dt = datetime(year, month, day, hour, minute, second, milli, utc)
        return dt, leap_second

    def utc_iso(self, places=0):
        """Return this UTC value as an ISO 8601 string or array of strings.

        For example: ``2014-01-18T01:35:38Z``

        """
        if places:
            power_of_ten = 10 ** places
            offset = _half_second / power_of_ten
            y, m, d, H, M, S = self._utc(offset)
            S, F = divmod(S, 1.0)
            format = '%%04d-%%02d-%%02dT%%02d:%%02d:%%02d.%%0%ddZ' % places
            args = (y, m, d, H, M, S, F * power_of_ten)
        else:
            format = '%04d-%02d-%02dT%02d:%02d:%02dZ'
            args = self._utc(_half_second)

        if self.shape:
            return [format % tup for tup in zip(*args)]
        else:
            return format % args

    def utc_jpl(self):
        """Return UTC in the format used by the JPL HORIZONS system.

        For example: ``A.D. 2014-Jan-18 01:35:37.5000 UT``

        """
        offset = _half_second / 1e4
        y, m, d, H, M, S = self._utc(offset)
        S, F = divmod(S, 1.0)
        # TODO: B.C.?
        format = 'A.D. %04d-%s-%02d %02d:%02d:%02d.%04d UT'
        args = (y, _months[m], d, H, M, S, F * 1e4)

        if self.shape:
            return [format % tup for tup in zip(*args)]
        else:
            return format % args

    def utc_strftime(self, format):
        """Format this UTC time according to a standard format string.

        This internally calls the Python ``strftime()`` routine from the
        Standard Library ``time()`` module, for which you can find a
        quick reference at http://strftime.org/.

        """
        tup = self._utc(_half_second)
        y, mon, d, h, m, s = tup
        zero = zeros_like(y)
        tup = (y.astype(int), mon.astype(int), d.astype(int),
               h.astype(int), m.astype(int), s.astype(int),
               zero, zero, zero)
        if self.shape:
            return [strftime(format, item) for item in zip(*tup)]
        else:
            return strftime(format, tup)

    def _utc(self, offset=0.0):
        """Return UTC as (year, month, day, hour, minute, second.fraction).

        The `offset` is added to the UTC time before it is split into
        its components.  This is useful if the user is going to round
        the result before displaying it.  If the result is going to be
        displayed as seconds, for example, set `offset` to half a second
        and then throw away the fraction; if the result is going to be
        displayed as minutes, set `offset` to thirty seconds and then
        throw away the seconds; and so forth.

        """
        tai = self.tai + offset
        leap_dates, leap_offsets = self.cache.run(usno_leapseconds)
        leap_reverse_dates = leap_dates + leap_offsets / DAY_S
        i = searchsorted(leap_reverse_dates, tai, 'right')
        j = tai - leap_offsets[i] / DAY_S
        whole, fraction = divmod(j + 0.5, 1.0)
        whole = whole.astype(int)
        y, mon, d = calendar_date(whole)
        h, hfrac = divmod(fraction * 24.0, 1.0)
        m, s = divmod(hfrac * 3600.0, 60.0)
        is_leap_second = j < leap_dates[i-1]
        s += is_leap_second
        self.utc = utc = (y, mon, d, h.astype(int), m.astype(int), s)
        return utc

    def __getattr__(self, name):

        # Cache of several expensive functions of time.

        if name == 'P':
            self.P = P = compute_precession(self.tdb)
            return P

        if name == 'PT':
            self.PT = PT = rollaxis(self.P, 1)
            return PT

        if name == 'N':
            self.N = N = compute_nutation(self)
            return N

        if name == 'NT':
            self.NT = NT = rollaxis(self.N, 1)
            return NT

        if name == 'M':
            self.M = M = einsum('ij...,jk...,kl...->il...', self.N, self.P, B)
            return M

        if name == 'MT':
            self.MT = MT = rollaxis(self.M, 1)
            return MT

        # Conversion between timescales.

        if name == 'tai':
            self.tai = tai = self.tt - tt_minus_tai
            return tai

        if name == 'utc':
            tai = self.tai
            leap_dates, leap_offsets = self.cache.run(usno_leapseconds)
            leap_reverse_dates = leap_dates + leap_offsets / DAY_S
            i = searchsorted(leap_reverse_dates, tai, 'right')
            j = tai - leap_offsets[i] / DAY_S
            y, mon, d, h = cal_date(j)
            h, hfrac = divmod(h, 1.0)
            m, s = divmod(hfrac * 3600.0, 60.0)
            tup = (y, mon, d, h, m, s)
            return array(tup) if self.shape else tup

        if name == 'tdb':
            tt = self.tt
            self.tdb = tdb = tt + tdb_minus_tt(tt) / DAY_S
            return tdb

        if name == 'ut1':
            self.ut1 = ut1 = self.tt - self.delta_t / DAY_S
            return ut1

        raise AttributeError('no such attribute %r' % name)

    def __eq__(self, other_jd):
        return self.tt == other_jd.tt


def now():
    """Return the current date and time as a `JulianDate` object.

    For the return value to be correct, your operating system time and
    timezone settings must be set so that the Python Standard Library
    constructor ``datetime.datetime.utcnow()`` returns a correct UTC
    date and time.

    """
    return JulianDate(utc=datetime.utcnow().replace(tzinfo=utc))

def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day integer."""
    janfeb = month < 3
    return (day
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075)

def julian_date(year, month=1, day=1, hour=0, minute=0, second=0.0):
    """Given a proleptic Gregorian calendar date, return a Julian date float."""
    return julian_day(year, month, day) - 0.5 + (
        second + minute * 60.0 + hour * 3600.0) / DAY_S

def cal_date(jd):
    """Given a Julian Date, return a Gregorian year, month, day, and hour."""
    jd = _to_array(jd) + 0.5

    hour = jd % 1.0 * 24.0
    k = jd.astype(int) + 68569
    n = 4 * k // 146097;

    k = k - (146097 * n + 3) // 4
    m = 4000 * (k + 1) // 1461001
    k = k - 1461 * m // 4 + 31
    month = 80 * k // 2447
    day = k - 2447 * month // 80
    k = month // 11

    month = month + 2 - 12 * k
    year = 100 * (n - 49) + m + k

    return year, month, day, hour

def calendar_date(jd_integer):
    """Convert Julian Day `jd_integer` into a Gregorian (year, month, day)."""

    k = jd_integer + 68569
    n = 4 * k // 146097

    k = k - (146097 * n + 3) // 4
    m = 4000 * (k + 1) // 1461001
    k = k - 1461 * m // 4 + 31
    month = 80 * k // 2447
    day = k - 2447 * month // 80
    k = month // 11

    month = month + 2 - 12 * k
    year = 100 * (n - 49) + m + k

    return year, month, day

def tdb_minus_tt(jd_tdb):
    """Computes TT corresponding to a TDB Julian date."""

    t = (jd_tdb - T0) / 36525.0

    # Expression given in USNO Circular 179, eq. 2.6.

    return (0.001657 * sin ( 628.3076 * t + 6.2401)
          + 0.000022 * sin ( 575.3385 * t + 4.2970)
          + 0.000014 * sin (1256.6152 * t + 6.1969)
          + 0.000005 * sin ( 606.9777 * t + 4.0212)
          + 0.000005 * sin (  52.9691 * t + 0.4444)
          + 0.000002 * sin (  21.3299 * t + 5.5431)
          + 0.000010 * t * sin ( 628.3076 * t + 4.2490))

def usno_leapseconds(cache):
    """Download the USNO table of leap seconds as a ``(2, N+1)`` NumPy array.

    The array has two rows ``[leap_dates leap_offsets]``.  The first row
    is used to find where a given date ``jd`` falls in the table::

        index = np.searchsorted(leap_dates, jd, 'right')

    This can return a value from ``0`` to ``N``, allowing the
    corresponding UTC offset to be fetched with::

        offset = leap_offsets[index]

    The offset is the number of seconds that must be added to a UTC time
    to build the corresponding TAI time.

    """
    with cache.open_url('http://maia.usno.navy.mil/ser7/leapsec.dat') as f:
        lines = f.readlines()

    linefields = [line.split() for line in lines]
    dates = [float(fields[4]) for fields in linefields]
    offsets = [float(fields[6]) for fields in linefields]

    dates.insert(0, float('-inf'))
    dates.append(float('inf'))

    offsets.insert(0, offsets[0])
    offsets.insert(1, offsets[0])

    return array([dates, offsets])

def _utc_datetime_to_tai(leap_dates, leap_offsets, dt):
    tup = dt.astimezone(utc).utctimetuple()
    year, month, day, hour, minute, second, wday, yday, dst = tup
    return _utc_to_tai(leap_dates, leap_offsets, year, month, day,
                       hour, minute, second + dt.microsecond / 1000000.00)

def _utc_date_to_tai(leap_dates, leap_offsets, d):
    return _utc_to_tai(leap_dates, leap_offsets, d.year, d.month, d.day)

def _utc_to_tai(leap_dates, leap_offsets, year, month=1, day=1,
                hour=0, minute=0, second=0.0):
    j = julian_day(year, month, day) - 0.5
    i = searchsorted(leap_dates, j, 'right')
    return j + (second + leap_offsets[i]
                + minute * 60.0
                + hour * 3600.0) / DAY_S
