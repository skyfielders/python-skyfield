from datetime import date, datetime, timedelta, tzinfo
from numpy import (array, concatenate, einsum, float_, interp, isnan, nan,
                   rollaxis, searchsorted, sin, where, zeros_like)
from time import strftime
from .constants import B1950, DAY_S, T0
from .earthlib import sidereal_time
from .framelib import ICRS_to_J2000 as B
from .functions import load_bundled_npy
from .nutationlib import compute_nutation, earth_tilt
from .precessionlib import compute_precession

try:
    from pytz import utc
except ImportError:

    # Lacking a full suite of timezones from pytz, we need to at least a
    # time zone object for UTC.

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

def _to_array(value):
    """When `value` is a plain Python sequence, return it as a NumPy array."""
    if hasattr(value, 'shape'):
        return value
    elif hasattr(value, '__len__'):
        return array(value)
    else:
        return float_(value)

tt_minus_tai = array(32.184 / DAY_S)

class Timescale(object):
    """The data necessary to express dates in different timescales.

    Whenever you want to express a date in Skyfield, you need a
    `Timescale` that can translate between several different systems for
    expressing time.  You will usually create a single `Timescale` at
    the beginning of your program, and use it every time you want to
    generate a specific `Time`:

    >>> from skyfield.api import load
    >>> ts = load.timescale()
    >>> t = ts.utc(1980, 3, 1, 9, 30)
    >>> t
    <Time tt=2444299.896426>

    Loading a timescale downloads tables from the United States Naval
    Observatory and the International Earth Rotation Service.  These
    files go out of date, and Skyfield will fetch updated copies once
    your copy of the files are too old.

    """
    _utcnow = datetime.utcnow

    def __init__(self, delta_t_recent, leap_dates, leap_offsets):
        self.delta_t_table = build_delta_t_table(delta_t_recent)
        self.leap_dates, self.leap_offsets = leap_dates, leap_offsets
        self.J2000 = Time(self, float_(T0))
        self.B1950 = Time(self, float_(B1950))

    def now(self):
        """Return the current date and time as a `Time` object.

        For the return value to be correct, your operating system time
        and timezone settings must be set so that the Python Standard
        Library constructor ``datetime.datetime.utcnow()`` returns a
        correct UTC date and time.

        """
        return self.utc(self._utcnow().replace(tzinfo=utc))

    def utc(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Return the Time corresponding to a specific moment in UTC.

        You can either specify the date as separate components, or
        provide a time zone aware Python datetime.  The following two
        calls are equivalent (the ``utc`` time zone object can be
        imported from the ``skyfield.api`` module, or from ``pytz`` if
        you have it)::

            ts.utc(2014, 1, 18, 1, 35, 37.5)
            ts.utc(datetime(2014, 1, 18, 1, 35, 37, 500000, tzinfo=utc))

        Note that only by passing the components separately can you
        specify a leap second, because a Python datetime will not allow
        the value 60 in its seconds field.

        """
        if isinstance(year, datetime):
            dt = year
            tai = _utc_datetime_to_tai(self.leap_dates, self.leap_offsets, dt)
        elif isinstance(year, date):
            d = year
            tai = _utc_date_to_tai(self.leap_dates, self.leap_offsets, d)
        elif hasattr(year, '__len__') and isinstance(year[0], datetime):
            # TODO: clean this up and better document the possibilities.
            list_of_datetimes = year
            tai = array([
                _utc_datetime_to_tai(self.leap_dates, self.leap_offsets, dt)
                for dt in list_of_datetimes])
        else:
            tai = _utc_to_tai(self.leap_dates, self.leap_offsets,
                              _to_array(year), _to_array(month),
                              _to_array(day), _to_array(hour),
                              _to_array(minute), _to_array(second))
        t = Time(self, tai + tt_minus_tai)
        t.tai = tai
        return t

    def tai(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
            jd=None):
        """Return the Time corresponding to a specific moment in TAI.

        You can specify International Atomic Time (TAI) by providing
        either a proleptic Gregorian calendar date or a raw Julian Date
        float.  The following two method calls are equivalent::

            timescale.tai(2014, 1, 18, 1, 35, 37.5)
            timescale.tai(jd=2456675.56640625)

        """
        if jd is not None:
            tai = jd
        else:
            tai = julian_date(year, month, day, hour, minute, second)
        tai = _to_array(tai)
        t = Time(self, tai + tt_minus_tai)
        t.tai = tai
        return t

    def tt(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
           jd=None):
        """Return the Time corresponding to a specific moment in TT.

        You can supply the Terrestrial Time (TT) by providing either a
        proleptic Gregorian calendar date or a raw Julian Date float.
        The following two method calls are equivalent::

            timescale.tt(2014, 1, 18, 1, 35, 37.5)
            timescale.tt(jd=2456675.56640625)

        """
        if jd is not None:
            tt = jd
        else:
            tt = julian_date(year, month, day, hour, minute, second)
        tt = _to_array(tt)
        return Time(self, tt)

    def tdb(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
            jd=None):
        """Return the Time corresponding to a specific moment in TDB.

        You can supply the Barycentric Dynamical Time (TDB) by providing
        either a proleptic Gregorian calendar date or a raw Julian Date
        float.  The following two method calls are equivalent::

            timescale.tdb(2014, 1, 18, 1, 35, 37.5)
            timescale.tdb(jd=2456675.56640625)

        """
        if jd is not None:
            tdb = jd
        else:
            tdb = julian_date(year, month, day, hour, minute, second)
        tdb = _to_array(tdb)
        tt = tdb - tdb_minus_tt(tdb) / DAY_S
        t = Time(self, tt)
        t.tdb = tdb
        return t

    def from_astropy(self, t):
        """Return a Skyfield time corresponding to the AstroPy time `t`."""
        return self.tt(jd=t.tt.jd)

class Time(object):
    """A single moment in history, or an array of several moments.

    You will typically not instantiate this class yourself, but will
    rely on a Skyfield ``Timescale`` object to build dates for you:

    >>> ts = load.timescale()
    >>> print(ts.utc(1980, 1, 1))
    <Time tt=2444239.500592>

    Times are represented internally by floating point Julian dates, but
    can be converted to other formats by using the many methods that
    time objects make available.

    """
    psi_correction = 0.0
    eps_correction = 0.0

    def __init__(self, ts, tt):
        self.tt = tt
        self.ts = ts
        self.shape = getattr(tt, 'shape', ())

    def __len__(self):
        return self.shape[0]

    def __repr__(self):

        if getattr(self.tt, 'shape', ()):
            rstr = '<Time {0} values from tt={1:.6f} to tt={2:.6f}>'
            return rstr.format(self.tt.size, self.tt.min(), self.tt.max())
        return '<Time tt={0:.6f}>'.format(self.tt)

    def __getitem__(self, index):
        # TODO: also copy cached matrices?
        t = Time(self.ts, self.tt[index])
        for name in 'tai', 'tdb', 'ut1', 'delta_t':
            value = getattr(self, name, None)
            if value is not None:
                if getattr(value, 'shape', None):
                    value = value[index]
                setattr(t, name, value)
        return t

    def astimezone(self, tz):
        """Convert to a Python ``datetime`` in a ``pytz`` provided timezone.

        Convert this time to a ``datetime`` in the timezone ``tz``,
        which should be one of the timezones provided by the third-party
        ``pytz`` package.  If this time is an array, then an array of
        datetimes is returned instead of a single value.

        """
        dt, leap_second = self.astimezone_and_leap_second(tz)
        return dt

    def astimezone_and_leap_second(self, tz):
        """Convert to a Python ``datetime`` and leap second in a timezone.

        Convert this time to a Python ``datetime`` and a leap second::

            dt, leap_second = t.astimezone_and_leap_second(tz)

        The argument ``tz`` should be a timezone from the third-party
        ``pytz`` package, which must be installed separately.  The date
        and time returned will be for that time zone.

        The leap second value is provided because a Python ``datetime``
        can only number seconds ``0`` through ``59``, but leap seconds
        have a designation of at least ``60``.  The leap second return
        value will normally be ``0``, but will instead be ``1`` if the
        date and time are a UTC leap second.  Add the leap second value
        to the ``second`` field of the ``datetime`` to learn the real
        name of the second.

        If this time is an array, then an array of ``datetime`` objects
        and an array of leap second integers is returned, instead of a
        single value each.

        """
        dt, leap_second = self.utc_datetime_and_leap_second()
        normalize = getattr(tz, 'normalize', None)
        if self.shape and normalize is not None:
            dt = array([normalize(d.astimezone(tz)) for d in dt])
        elif self.shape:
            dt = array([d.astimezone(tz) for d in dt])
        elif normalize is not None:
            dt = normalize(dt.astimezone(tz))
        else:
            dt = dt.astimezone(tz)
        return dt, leap_second

    def toordinal(self):
        """Return the proleptic Gregorian ordinal of the UTC date.

        This method makes Skyfield `Time` objects compatible with Python
        `datetime`_ objects, which also provide a ``toordinal()``
        method.  Thanks to this method, a `Time` can often be used
        directly as a coordinate for a plot.

        """
        return self._utc_float() - 1721424.5

    def utc_datetime(self):
        """Convert to a Python ``datetime`` in UTC.

        If the third-party `pytz`_ package is available, then its
        ``utc`` timezone will be used as the timezone of the returned
        `datetime`_.  Otherwise, an equivalent Skyfield ``utc`` timezone
        object is used.  If this time is an array, then a sequence of
        datetimes is returned instead of a single value.

        """
        dt, leap_second = self.utc_datetime_and_leap_second()
        return dt

    def utc_datetime_and_leap_second(self):
        """Convert to a Python ``datetime`` in UTC, plus a leap second value.

        Convert this time to a `datetime`_ object and a leap second::

            dt, leap_second = t.utc_datetime_and_leap_second()

        If the third-party `pytz`_ package is available, then its
        ``utc`` timezone will be used as the timezone of the return
        value.  Otherwise, Skyfield uses its own ``utc`` timezone.

        The leap second value is provided because a Python ``datetime``
        can only number seconds ``0`` through ``59``, but leap seconds
        have a designation of at least ``60``.  The leap second return
        value will normally be ``0``, but will instead be ``1`` if the
        date and time are a UTC leap second.  Add the leap second value
        to the ``second`` field of the ``datetime`` to learn the real
        name of the second.

        If this time is an array, then an array of ``datetime`` objects
        and an array of leap second integers is returned, instead of a
        single value each.

        """
        year, month, day, hour, minute, second = self._utc_tuple(
            _half_millisecond)
        second, fraction = divmod(second, 1.0)
        second = second.astype(int)
        leap_second = second // 60
        second -= leap_second
        milli = (fraction * 1000).astype(int) * 1000
        if self.shape:
            utcs = [utc] * self.shape[0]
            argsets = zip(year, month, day, hour, minute, second, milli, utcs)
            dt = array([datetime(*args) for args in argsets])
        else:
            dt = datetime(year, month, day, hour, minute, second, milli, utc)
        return dt, leap_second

    def utc_iso(self, places=0):
        """Convert to an ISO 8601 string like ``2014-01-18T01:35:38Z`` in UTC.

        If this time is an array of dates, then a sequence of strings is
        returned instead of a single string.

        """
        if places:
            power_of_ten = 10 ** places
            offset = _half_second / power_of_ten
            year, month, day, hour, minute, second = self._utc_tuple(offset)
            second, fraction = divmod(second, 1.0)
            fraction *= power_of_ten
            format = '%%04d-%%02d-%%02dT%%02d:%%02d:%%02d.%%0%ddZ' % places
            args = (year, month, day, hour, minute, second, fraction)
        else:
            format = '%04d-%02d-%02dT%02d:%02d:%02dZ'
            args = self._utc_tuple(_half_second)

        if self.shape:
            return [format % tup for tup in zip(*args)]
        else:
            return format % args

    def utc_jpl(self):
        """Convert to an ``A.D. 2014-Jan-18 01:35:37.5000 UT`` string.

        Returns a string for this date and time in UTC, in the format
        used by the JPL HORIZONS system.  If this time is an array of
        dates, then a sequence of strings is returned instead of a
        single string.

        """
        offset = _half_second / 1e4
        year, month, day, hour, minute, second = self._utc_tuple(offset)
        second, fraction = divmod(second, 1.0)
        fraction *= 1e4
        bc = year < 1
        year = abs(year - bc)
        era = where(bc, 'B.C.', 'A.D.')
        format = '%s %04d-%s-%02d %02d:%02d:%02d.%04d UT'
        args = (era, year, _months[month], day, hour, minute, second, fraction)

        if self.shape:
            return [format % tup for tup in zip(*args)]
        else:
            return format % args

    def utc_strftime(self, format):
        """Format the UTC time using a Python date formatting string.

        This internally calls the Python ``strftime()`` routine from the
        Standard Library ``time()`` module, for which you can find a
        quick reference at ``http://strftime.org/``.  If this object is
        an array of times, then a sequence of strings is returned
        instead of a single string.

        """
        tup = self._utc_tuple(_half_second)
        year, month, day, hour, minute, second = tup
        second = second.astype(int)
        zero = zeros_like(year)
        tup = (year, month, day, hour, minute, second, zero, zero, zero)
        if self.shape:
            return [strftime(format, item) for item in zip(*tup)]
        else:
            return strftime(format, tup)

    def _utc_tuple(self, offset=0.0):
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
        leap_dates = self.ts.leap_dates
        leap_offsets = self.ts.leap_offsets
        leap_reverse_dates = leap_dates + leap_offsets / DAY_S
        i = searchsorted(leap_reverse_dates, tai, 'right')
        j = tai - leap_offsets[i] / DAY_S
        whole, fraction = divmod(j + 0.5, 1.0)
        whole = whole.astype(int)
        year, month, day = calendar_date(whole)
        hour, hfrac = divmod(fraction * 24.0, 1.0)
        minute, second = divmod(hfrac * 3600.0, 60.0)
        is_leap_second = j < leap_dates[i-1]
        second += is_leap_second
        return year, month, day, hour.astype(int), minute.astype(int), second

    def _utc_float(self):
        """Return UTC as a floating point Julian date."""
        tai = self.tai
        leap_dates = self.ts.leap_dates
        leap_offsets = self.ts.leap_offsets
        leap_reverse_dates = leap_dates + leap_offsets / DAY_S
        i = searchsorted(leap_reverse_dates, tai, 'right')
        return tai - leap_offsets[i] / DAY_S

    def tai_calendar(self):
        """Return TAI as a tuple (year, month, day, hour, minute, second)."""
        return calendar_tuple(self.tai)

    def tt_calendar(self):
        """Return TT as a tuple (year, month, day, hour, minute, second)."""
        return calendar_tuple(self.tt)

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
            utc = self._utc_tuple()
            utc = array(utc) if self.shape else utc
            self.utc = utc = utc
            return utc

        if name == 'tdb':
            tt = self.tt
            self.tdb = tdb = tt + tdb_minus_tt(tt) / DAY_S
            return tdb

        if name == 'ut1':
            self.ut1 = ut1 = self.tt - self.delta_t / DAY_S
            return ut1

        if name == 'delta_t':
            table = self.ts.delta_t_table
            self.delta_t = delta_t = interpolate_delta_t(table, self.tt)
            return delta_t

        if name == 'gmst':
            self.gmst = gmst = sidereal_time(self)
            return gmst

        if name == 'gast':
            self.gast = gast = self.gmst + earth_tilt(self)[2] / 3600.0
            return gast

        raise AttributeError('no such attribute %r' % name)

    def __eq__(self, other_time):
        if not isinstance(other_time, Time):
            return NotImplemented
        return self.tt == other_time.tt

def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
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

def calendar_tuple(jd_float, offset=0.0):
    """Return a (year, month, day, hour, minute, second.fraction) tuple.

    The `offset` is added to the time before it is split into its
    components.  This is useful if the user is going to round the
    result before displaying it.  If the result is going to be
    displayed as seconds, for example, set `offset` to half a second
    and then throw away the fraction; if the result is going to be
    displayed as minutes, set `offset` to thirty seconds and then
    throw away the seconds; and so forth.

    """
    whole, fraction = divmod(jd_float + 0.5, 1.0)
    whole = whole.astype(whole)
    year, month, day = calendar_date(whole)
    hour, hfrac = divmod(fraction * 24.0, 1.0)
    minute, second = divmod(hfrac * 3600.0, 60.0)
    return year, month, day, hour.astype(int), minute.astype(int), second

def tdb_minus_tt(jd_tdb):
    """Computes how far TDB is in advance of TT, given TDB.

    Given that the two time scales never diverge by more than 2ms, TT
    can also be given as the argument to perform the conversion in the
    other direction.

    """
    t = (jd_tdb - T0) / 36525.0

    # USNO Circular 179, eq. 2.6.
    return (0.001657 * sin ( 628.3076 * t + 6.2401)
          + 0.000022 * sin ( 575.3385 * t + 4.2970)
          + 0.000014 * sin (1256.6152 * t + 6.1969)
          + 0.000005 * sin ( 606.9777 * t + 4.0212)
          + 0.000005 * sin (  52.9691 * t + 0.4444)
          + 0.000002 * sin (  21.3299 * t + 5.5431)
          + 0.000010 * t * sin ( 628.3076 * t + 4.2490))

def interpolate_delta_t(delta_t_table, tt):
    """Return interpolated Delta T values for the times in `tt`.

    The 2xN table should provide TT values as element 0 and
    corresponding Delta T values for element 1.  For times outside the
    range of the table, a long-term formula is used instead.

    """
    tt_array, delta_t_array = delta_t_table
    delta_t = _to_array(interp(tt, tt_array, delta_t_array, nan, nan))
    missing = isnan(delta_t)

    if missing.any():
        # Test if we are dealing with an array and proceed appropriately
        if missing.shape:
            tt = tt[missing]
            delta_t[missing] = delta_t_formula_morrison_and_stephenson_2004(tt)
        else:
            delta_t = delta_t_formula_morrison_and_stephenson_2004(tt)
    return delta_t

def delta_t_formula_morrison_and_stephenson_2004(tt):
    """Delta T formula from Morrison and Stephenson, 2004.

    This parabola can be used to estimate the value of Delta T for dates
    in the far past or future, for which more specific estimates are not
    available.

    """
    t = (tt - 2385800.5) / 36525.0  # centuries before or after 1820
    return 32.0 * t * t - 20.0

def build_delta_t_table(delta_t_recent):
    """Build a table for interpolating Delta T.

    Given a 2xN array of recent Delta T values, whose element 0 is a
    sorted array of TT Julian dates and element 1 is Delta T values,
    this routine returns a more complete table by prepending two
    built-in data sources that ship with Skyfield as pre-built arrays:

    * The historical values from Morrison and Stephenson (2004) which
      the http://eclipse.gsfc.nasa.gov/SEcat5/deltat.html NASA web page
      presents in an HTML table.

    * The United States Naval Observatory ``historic_deltat.data``
      values for Delta T over the years 1657 through 1984.

    """
    ancient = load_bundled_npy('morrison_stephenson_deltat')
    historic = load_bundled_npy('historic_deltat')

    # Prefer USNO over Morrison and Stephenson where they overlap.
    historic_start_time = historic[0,0]
    i = searchsorted(ancient[0], historic_start_time)
    bundled = concatenate([ancient[:,:i], historic], axis=1)

    # Let recent data replace everything else.
    recent_start_time = delta_t_recent[0,0]
    i = searchsorted(bundled[0], recent_start_time)
    row = ((0,),(0,))
    table = concatenate([row, bundled[:,:i], delta_t_recent, row], axis=1)

    # Create initial and final point to provide continuity with formula.
    century = 36524.0
    start = table[0,1] - century
    table[:,0] = start, delta_t_formula_morrison_and_stephenson_2004(start)
    end = table[0,-2] + century
    table[:,-1] = end, delta_t_formula_morrison_and_stephenson_2004(end)
    return table

def _utc_datetime_to_tai(leap_dates, leap_offsets, dt):
    try:
        utc_datetime = dt.astimezone(utc)
    except ValueError:
        raise ValueError(_naive_complaint)
    tup = utc_datetime.utctimetuple()
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

_JulianDate_deprecation_message = """Skyfield no longer supports direct\
 instantiation of JulianDate objects (which are now called Time objects)

If you need to quickly get an old Skyfield script working again, simply
downgrade to Skyfield version 0.6.1 using a command like:

        pip install skyfield==0.6.1

Skyfield used to let you build JulianDate objects directly:

        t = JulianDate(utc=(1980, 4, 20))   # the old way

But this forced Skyfield to maintain secret global copies of several
time scale data files, that need to be downloaded and kept up to date
for Time objects to work.  Skyfield now makes this collection of data
files explicit, and calls the bundle of files a "Timescale" object.  You
can create one with the "load.timescale()" method and then build times
using its methods, which let you either specify a calendar date or else
supply a raw Julian date value with the "jd" keyword:

        from skyfield.api import load
        ts = load.timescale()
        t = ts.utc(1980, 4, 20)       # the new way

        t = ts.tt(jd=2444349.500592)  # jd is also supported for tai, tt, tdb

See http://rhodesmill.org/skyfield/time.html for more details."""


_naive_complaint = """cannot interpret a datetime that lacks a timezone

You must either specify that your datetime is in UTC:

    from skyfield.api import utc
    d = datetime(..., tzinfo=utc)  # to build a new datetime
    d = d.replace(tzinfo=utc)      # to fix an existing datetime

Or install the third-party `pytz` library and use any of its timezones:

    from pytz import timezone
    eastern = timezone('US/Eastern')
    d = eastern.localize(datetime(2014, 1, 16, 1, 32, 9))"""
