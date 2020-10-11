# -*- coding: utf-8 -*-
import datetime as dt_module
import re
from collections import namedtuple
from datetime import date, datetime
from numpy import (array, concatenate, cos, float_, interp, isnan, nan,
                   ndarray, pi, rollaxis, searchsorted, sin, where, zeros_like)
from time import strftime, struct_time
from .constants import ASEC2RAD, B1950, DAY_S, T0, tau
from .descriptorlib import reify
from .earthlib import sidereal_time, earth_rotation_angle
from .framelib import ICRS_to_J2000 as B
from .functions import (mxm, mxmxm, load_bundled_npy, rot_z,
                        _to_array, _reconcile)
from .nutationlib import (
    build_nutation_matrix, equation_of_the_equinoxes_complimentary_terms,
    iau2000a_radians, mean_obliquity,
)
from .precessionlib import compute_precession

GREGORIAN_START = 2299161
GREGORIAN_START_ENGLAND = 2361222

CalendarTuple = namedtuple('CalendarTuple', 'year month day hour minute second')

class CalendarArray(ndarray):
    @property
    def year(self): return self[0]
    @property
    def month(self): return self[1]
    @property
    def day(self): return self[2]
    @property
    def hour(self): return self[3]
    @property
    def minute(self): return self[4]
    @property
    def second(self): return self[5]

if hasattr(dt_module, 'timezone'):
    utc = dt_module.timezone.utc
else:
    class UTC(dt_module.tzinfo):
        'UTC'
        zero = dt_module.timedelta(0)
        def utcoffset(self, dt):
            return self.zero
        def tzname(self, dt):
            return 'UTC'
        def dst(self, dt):
            return self.zero

    utc = UTC()

# Much of the following code is adapted from the USNO's "novas.c".

_time_zero = dt_module.time(tzinfo=utc)
_half_minute = 30.0 / DAY_S
_half_second = 0.5 / DAY_S
_half_microsecond = 0.5e-6 / DAY_S
_months = array(['Month zero', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

tt_minus_tai = array(32.184 / DAY_S)

class Timescale(object):
    """The data necessary to express dates in different timescales.

    A `Timescale` loads three data files.  They supply the value of ∆T
    over time and the schedule of UTC leap seconds, which Skyfield uses
    to translate between different time scales.  Most programs create a
    single `Timescale` which they use to build their `Time` objects:

    >>> from skyfield.api import load
    >>> ts = load.timescale()
    >>> t = ts.utc(1980, 3, 1, 9, 30)
    >>> t
    <Time tt=2444299.896425741>

    See :ref:`downloading-timescale-files` if you are interested in
    checking how recent the data is in the files loaded by the
    timescale.

    """
    _utcnow = datetime.utcnow

    def __init__(self, delta_t_recent, leap_dates, leap_offsets):
        self.delta_t_table = build_delta_t_table(delta_t_recent)
        self.leap_dates, self.leap_offsets = leap_dates, leap_offsets
        self._leap_reverse_dates = leap_dates + leap_offsets / DAY_S
        self.J2000 = Time(self, float_(T0))
        self.B1950 = Time(self, float_(B1950))
        self.julian_calendar_cutoff = None

    def now(self):
        """Return the current date and time as a `Time` object.

        For the return value to be correct, your operating system time
        and timezone settings must be set so that the Python Standard
        Library constructor ``datetime.datetime.utcnow()`` returns a
        correct UTC date and time.

        """
        return self.from_datetime(self._utcnow().replace(tzinfo=utc))

    def from_datetime(self, datetime):
        """Return a `Time` for a Python ``datetime``.

        The ``datetime`` must be “timezone-aware”: it must have a time
        zone object as its ``tzinfo`` attribute instead of ``None``.

        """
        return self._utc(_datetime_to_utc_tuple(datetime))

    def from_datetimes(self, datetime_list):
        """Return a `Time` for a list of Python ``datetime``\ s.

        The ``datetime`` objects must each be “timezone-aware”: they
        must each have a time zone object as their ``tzinfo`` attribute
        instead of ``None``.

        """
        tuples = (_datetime_to_utc_tuple(d) for d in datetime_list)
        return self._utc(array(value) for value in zip(*tuples))

    def utc(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Build a `Time` from a UTC `calendar date`."""
        # TODO: someday deprecate passing datetime objects here, as
        # there are now separate constructors for them.
        if isinstance(year, datetime):
            return self.from_datetime(year)
        if isinstance(year, date):
            return self.from_datetime(datetime.combine(year, _time_zero))
        if hasattr(year, '__len__') and isinstance(year[0], datetime):
            return self.from_datetimes(year)

        tup = year, month, day, hour, minute, second
        return self._utc(tup)

    def _utc(self, tup):
        year, month, day, hour, minute, second = tup
        whole, fraction = self._jd(year, month, day, hour, minute, 0.0)
        i = searchsorted(self.leap_dates, whole + fraction, 'right')
        fraction += (self.leap_offsets[i] + second) / DAY_S
        whole, fraction = _reconcile(whole, fraction)  # second could be array
        t = Time(self, whole, fraction + tt_minus_tai)
        t.tai_fraction = fraction
        return t

    def _jd(self, year, month, day, hour, minute, second):
        a = _to_array
        cutoff = self.julian_calendar_cutoff
        whole = julian_day(a(year), a(month), a(day), cutoff) - 0.5
        fraction = (a(second) + a(minute) * 60.0 + a(hour) * 3600.0) / DAY_S
        return _reconcile(whole, fraction)

    def _cal(self, whole, fraction):
        return calendar_tuple(whole, fraction, self.julian_calendar_cutoff)

    def _strftime(self, format, jd, fraction, seconds_bump=None):
        # Python forces an unhappy choice upon us: either use the faster
        # time.strftime() and lose support for '%f', or use the slower
        # datetime.strftime() and crash if years are negative.  We take the
        # first option, but then patch '%f' support back in by secretly
        # passing the microseconds string as the time zone name.  After all,
        # the routines supported by this function never use time zones.
        # What could go wrong?

        ms = _format_uses_milliseconds(format)

        if ms:
            fraction = fraction + 1e-16  # encourage .0 to not turn into .999999
        elif _format_uses_seconds(format):
            fraction = fraction + _half_second
        elif _format_uses_minutes(format):
            fraction = fraction + _half_minute

        year, month, day, hour, minute, second = self._cal(jd, fraction)
        z = year * 0

        # TODO: will this always produce the same whole number that
        # calendar_tuple() produces internally?  Or should we make a private
        # version of calendar_tuple() that returns it to us for use here?
        weekday = (fraction + 0.5 + _to_array(jd)).astype(int) % 7

        if _format_uses_day_of_year(format):
            start_of_year = julian_day(year, 1, 1, self.julian_calendar_cutoff)
            yday = (jd + fraction + 1.5 - start_of_year).astype(int)
        else:
            yday = z

        if ms:
            format = format[:ms.start()] + '%Z' + format[ms.end():]
            second = (second * 1e6).astype(int)
            second, usec = divmod(second, 1000000)
            if seconds_bump is not None:
                second += seconds_bump
            if getattr(jd, 'ndim', 0):
                u = ['%06d' % u for u in usec]
                tup = (year, month, day, hour, minute, second,
                       weekday, yday, z, u)
                return [strftime(format, struct_time(t)) for t in zip(*tup)]
            u = '%06d' % usec
            tup = year, month, day, hour, minute, second, weekday, yday, z, u
            return strftime(format, struct_time(tup))
        else:
            second = second.astype(int)
            if seconds_bump is not None:
                second += seconds_bump
            tup = year, month, day, hour, minute, second, weekday, yday, z
            if getattr(jd, 'ndim', 0):
                return [strftime(format, item) for item in zip(*tup)]
            return strftime(format, tup)

    def tai(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
            jd=None):
        """Build a `Time` from an International Atomic Time `calendar date`."""
        if jd is not None:
            return self.tai_jd(jd)  # deprecate someday
        whole, fraction = self._jd(year, month, day, hour, minute, second)
        t = Time(self, whole, fraction + tt_minus_tai)
        t.tai_fraction = fraction
        return t

    def tai_jd(self, jd, fraction=None):
        """Build a `Time` from an International Atomic Time Julian date."""
        jd, fraction = _normalize_jd_and_fraction(jd, fraction)
        t = Time(self, jd, fraction + tt_minus_tai)
        t.tai_fraction = fraction
        return t

    def tt(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
           jd=None):
        """Build a `Time` from a Terrestrial Time `calendar date`."""
        if jd is not None:
            return self.tt_jd(jd)  # deprecate someday
        whole, fraction = self._jd(year, month, day, hour, minute, second)
        return Time(self, whole, fraction)

    def tt_jd(self, jd, fraction=None):
        """Build a `Time` from a Terrestrial Time Julian date."""
        jd, fraction = _normalize_jd_and_fraction(jd, fraction)
        return Time(self, jd, fraction)

    def J(self, year):
        """Build a `Time` from a Terrestrial Time Julian year or array.

        Julian years are convenient uniform periods of exactly 365.25
        days of Terrestrial Time, centered on 2000 January 1 12h TT =
        Julian year 2000.0.

        """
        tt = _to_array(year) * 365.25 + 1721045.0
        return Time(self, tt, 0.0)

    def tdb(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
            jd=None):
        """Build a `Time` from a Barycentric Dynamical Time `calendar date`."""
        if jd is not None:
            return self.tdb_jd(jd)  # deprecate someday
        whole, fraction = self._jd(year, month, day, hour, minute, second)
        jd = whole + fraction  # TODO: why do tests break if we pass separately
        return Time(self, jd, - tdb_minus_tt(jd) / DAY_S)

    def tdb_jd(self, jd, fraction=None):
        """Build `Time` from a Barycentric Dynamical Time Julian date."""
        jd, fraction = _normalize_jd_and_fraction(jd, fraction)
        t = Time(self, jd, fraction - tdb_minus_tt(jd, fraction) / DAY_S)
        t.tdb_fraction = fraction
        return t

    def ut1(self, year=None, month=1, day=1, hour=0, minute=0, second=0.0,
            jd=None):
        """Build a `Time` from a UT1 Universal Time `calendar date`."""
        if jd is None:  # TODO: deprecate the jd parameter to this method
            whole, fraction = self._jd(year, month, day, hour, minute, second)
            jd = whole + fraction  # TODO: can we pass high precision on?
        return self.ut1_jd(jd)

    def ut1_jd(self, jd):
        """Build a `Time` from a UT1 Universal Time Julian date."""
        ut1 = _to_array(jd)

        # Estimate TT = UT1, to get a rough Delta T estimate.
        tt_approx = ut1
        delta_t_approx = interpolate_delta_t(self.delta_t_table, tt_approx)

        # Use the rough Delta T to make a much better estimate of TT,
        # then generate an even better Delta T.
        tt_approx = ut1 + delta_t_approx / DAY_S
        delta_t_approx = interpolate_delta_t(self.delta_t_table, tt_approx)

        # We can now estimate TT with an error of < 1e-9 seconds within
        # 10 centuries of either side of the present; for details, see:
        # https://github.com/skyfielders/astronomy-notebooks
        # and look for the notebook "error-in-timescale-ut1.ipynb".
        delta_t_approx /= DAY_S
        t = Time(self, ut1, delta_t_approx)
        t.ut1_fraction = 0.0 * ut1
        return t

    def from_astropy(self, t):
        """Build a Skyfield `Time` from an AstroPy time object."""
        return self.tt(jd=t.tt.jd)

class Time(object):
    """A single moment in history, or an array of several moments.

    Skyfield programs don’t usually instantiate this class directly, but
    instead build time objects using one of the timescale methods listed
    at `timescale-summary`.  If you do attempt the low-level operation
    of building a time object yourself, either leave ``tt_fraction`` at
    its default value of ``None`` — in which case Skyfield will assume
    the fraction is zero — or provide a ``tt_fraction`` array that has
    exactly the same dimensions as your ``tt`` array.

    """
    psi_correction = 0.0  # TODO: temporarily unsupported
    eps_correction = 0.0  # TODO: temporarily unsupported

    def __init__(self, ts, tt, tt_fraction=None):
        if tt_fraction is None:
            tt_fraction = zeros_like(tt)
        self.ts = ts
        self.whole = tt
        self.tt_fraction = tt_fraction
        self.shape = getattr(tt, 'shape', ())

    def __len__(self):
        return self.shape[0]

    def __repr__(self):
        size = getattr(self.tt, 'size', -1)
        if size > 3:
            rstr = '<Time tt=[{0} ... {1}] len={2}>'
            return rstr.format(self.tt[0], self.tt[-1], size)
        else:
            return ('<Time tt={0}>'.format(self.tt)
                    .replace('[ ', '[').replace('  ', ' '))

    def __getitem__(self, index):
        # TODO: also copy cached matrices?
        # TODO: raise non-IndexError exception if this Time is not an array;
        # otherwise, a `for` loop over it will not raise an error.
        t = Time(self.ts, self.whole[index], self.tt_fraction[index])
        for name in 'tai_fraction', 'tdb_fraction', 'ut1_fraction':
            # TODO: drat, I suspect this forces the creation of all of
            # the fractions whether we need them or not.
            value = getattr(self, name, None)
            if value is not None:
                setattr(t, name, value[index])
        return t

    def astimezone(self, tz):
        """Convert to a Python ``datetime`` in a particular timezone ``tz``.

        If this time is an array, then an array of datetimes is returned
        instead of a single value.

        """
        dt, leap_second = self.astimezone_and_leap_second(tz)
        return dt

    def astimezone_and_leap_second(self, tz):
        """Convert to a Python ``datetime`` and leap second in a timezone.

        Convert this time to a Python ``datetime`` and a leap second::

            dt, leap_second = t.astimezone_and_leap_second(tz)

        The argument ``tz`` should be a ``datetime`` compatible
        timezone.

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
        whole, fraction, is_leap_second = self._utc_float(0.0)
        return whole - 1721424.5 + fraction

    def utc_datetime(self):
        """Convert to a Python ``datetime`` in UTC.

        If this time is an array, then a list of datetimes is returned
        instead of a single value.

        """
        dt, leap_second = self.utc_datetime_and_leap_second()
        return dt

    def utc_datetime_and_leap_second(self):
        """Convert to a Python ``datetime`` in UTC, plus a leap second value.

        Convert this time to a `datetime`_ object and a leap second::

            dt, leap_second = t.utc_datetime_and_leap_second()

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
            _half_microsecond)
        micro = (second * 1e6).astype(int)
        second, micro = divmod(micro, 1000000)
        leap_second = second // 60
        second -= leap_second  # fit within limited bounds of Python datetime
        if self.shape:
            zone = [utc] * self.shape[0]
            argsets = zip(year, month, day, hour, minute, second, micro, zone)
            dt = array([datetime(*args) for args in argsets])
        else:
            dt = datetime(year, month, day, hour, minute, second, micro, utc)
        return dt, leap_second

    def utc_iso(self, delimiter='T', places=0):
        """Convert to an ISO 8601 string like ``2014-01-18T01:35:38Z`` in UTC.

        If this time is an array of dates, then a sequence of strings is
        returned instead of a single string.

        """
        # "places" used to be the 1st argument, so continue to allow an
        # integer in that spot.  TODO: deprecate this in Skyfield 2.0
        # and remove it in 3.0.
        if isinstance(delimiter, int):
            places = delimiter
            delimiter = 'T'

        if places:
            power_of_ten = 10 ** places
            offset = _half_second / power_of_ten
            year, month, day, hour, minute, second = self._utc_tuple(offset)
            second, fraction = divmod(second, 1.0)
            fraction *= power_of_ten
            format = '%04d-%02d-%02d{0}%02d:%02d:%02d.%0{1}dZ'.format(
                delimiter, places)
            args = (year, month, day, hour, minute, second, fraction)
        else:
            format = '%04d-%02d-%02d{0}%02d:%02d:%02dZ'.format(delimiter)
            args = self._utc_tuple(_half_second)

        if self.shape:
            return [format % tup for tup in zip(*args)]
        else:
            return format % args

    def utc_jpl(self):
        """Convert to a string like ``A.D. 2014-Jan-18 01:35:37.5000 UT``.

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

    def utc_strftime(self, format='%Y-%m-%d %H:%M:%S UTC'):
        """Format the UTC time using a Python datetime formatting string.

        This calls Python’s ``time.strftime()`` to format the date and
        time.  A single string is returned or else a whole array of
        strings, depending on whether this time object is an array.
        The most commonly used formats are:

        * ``%Y`` four-digit year, ``%y`` two-digit year
        * ``%m`` month number, ``%B`` name, ``%b`` abbreviation
        * ``%d`` day of month
        * ``%H`` hour
        * ``%M`` minute
        * ``%S`` second
        * ``%A`` day of week, ``%a`` its abbreviation

        You can find the full list, along with options that control
        field widths and leading zeros, at:

        https://docs.python.org/3/library/time.html#time.strftime

        If the smallest time unit in your format is minutes or seconds,
        then the time is rounded to the nearest minute or second.
        Otherwise the value is truncated rather than rounded.

        """
        whole, fraction, is_leap_second = self._utc_float(0.0)
        return self.ts._strftime(format, self.whole, fraction, is_leap_second)

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
        jd, fraction, is_leap_second = self._utc_float(offset)
        year, month, day, hour, minute, second = self.ts._cal(jd, fraction)
        second += is_leap_second
        return year, month, day, hour.astype(int), minute.astype(int), second

    def _utc_float(self, offset):
        ts = self.ts
        i = searchsorted(ts._leap_reverse_dates, self.tai + offset, 'right')
        whole = self.whole
        fraction = offset - ts.leap_offsets[i] / DAY_S + self.tai_fraction
        is_leap_second = (whole + fraction) < ts.leap_dates[i-1]
        return whole, fraction, is_leap_second

    # Calendar tuples.

    def tai_calendar(self):
        """TAI as a (year, month, day, hour, minute, second) `calendar date`."""
        return self.ts._cal(self.whole, self.tai_fraction)

    def tt_calendar(self):
        """TT as a (year, month, day, hour, minute, second) `calendar date`."""
        return self.ts._cal(self.whole, self.tt_fraction)

    def tdb_calendar(self):
        """TDB as a (year, month, day, hour, minute, second) `calendar date`."""
        return self.ts._cal(self.whole, self.tdb_fraction)

    def ut1_calendar(self):
        """UT1 as a (year, month, day, hour, minute, second) `calendar date`."""
        return self.ts._cal(self.whole, self.ut1_fraction)

    # Date formatting.

    def tai_strftime(self, format='%Y-%m-%d %H:%M:%S TAI'):
        """Format TAI with a datetime strftime() format string."""
        return self.ts._strftime(format, self.whole, self.tai_fraction)

    def tt_strftime(self, format='%Y-%m-%d %H:%M:%S TT'):
        """Format TT with a datetime strftime() format string."""
        return self.ts._strftime(format, self.whole, self.tt_fraction)

    def tdb_strftime(self, format='%Y-%m-%d %H:%M:%S TDB'):
        """Format TDB with a datetime strftime() format string."""
        return self.ts._strftime(format, self.whole, self.tdb_fraction)

    def ut1_strftime(self, format='%Y-%m-%d %H:%M:%S UT1'):
        """Format UT1 with a datetime strftime() format string."""
        return self.ts._strftime(format, self.whole, self.ut1_fraction)

    # Convenient caching of several expensive functions of time.

    @reify
    def M(self):
        """3×3 rotation matrix: ICRS → equinox of this date."""

        # Compute N and P instead of asking for self.N and self.P to
        # avoid keeping copies of them, since Skyfield itself never uses
        # them again once M has been computed.  But we do check in case
        # a user has forced them to be built already.

        d = self.__dict__

        P = d.get('P')
        if P is None:
            P = self.precession_matrix()

        N = d.get('N')
        if N is None:
            N = self.nutation_matrix()

        return mxmxm(N, P, B)

    @reify
    def MT(self):
        """3×3 rotation matrix: equinox of this date → ICRS."""
        return rollaxis(self.M, 1)

    @reify
    def C(self):
        # Calculate the Equation of Origins in cycles
        eq_origins = (earth_rotation_angle(self.ut1) - self.gast / 24.0)
        R = rot_z(2 * pi * eq_origins)
        return mxm(R, self.M)

    @reify
    def CT(self):
        return rollaxis(self.C, 1)

    @reify
    def _nutation_angles_radians(self):
        # TODO: add psi and eps corrections support back in here, rather
        # than at points of use.
        return iau2000a_radians(self)

    def _nutation_angles(self, angles):
        # Sample code shared with early adopters suggested that setting
        # this attribute manually could avoid the expense of IAU 2000A,
        # so this setter continues to support the pattern.

        d_psi, d_eps = angles
        self._nutation_angles_radians = (
            d_psi / 1e7 * ASEC2RAD,
            d_eps / 1e7 * ASEC2RAD,
        )

    _nutation_angles = property(None, _nutation_angles)

    @reify
    def _mean_obliquity_radians(self):
        # Cached because it is used to compute both gast and N.
        return mean_obliquity(self.tdb) * ASEC2RAD

    # Conversion between timescales.

    @reify
    def J(self):
        """Return a floating point Julian year or array of years for this date.

        Julian years are convenient uniform periods of exactly 365.25
        days of Terrestrial Time, centered on 2000 January 1 12h TT =
        Julian year 2000.0.

        """
        return (self.whole - 1721045.0 + self.tt_fraction) / 365.25

    @reify
    def utc(self):
        utc = self._utc_tuple()
        return (array(utc).view(CalendarArray) if self.shape
                else CalendarTuple(*utc))

    @reify
    def tai_fraction(self):
        return self.tt_fraction - tt_minus_tai

    @reify
    def tdb_fraction(self):
        fr = self.tt_fraction
        return fr + tdb_minus_tt(self.whole, fr) / DAY_S

    @reify
    def ut1_fraction(self):
        # Calling "self.delta_t" would cache a useless intermediate value, so:
        table = self.ts.delta_t_table
        delta_t = interpolate_delta_t(table, self.tt)
        return self.tt_fraction - delta_t / DAY_S

    @reify
    def delta_t(self):
        table = self.ts.delta_t_table
        return interpolate_delta_t(table, self.tt)

    @reify
    def dut1(self):
        ts = self.ts
        i = searchsorted(ts._leap_reverse_dates, self.tai, 'right')
        return 32.184 + ts.leap_offsets[i] - self.delta_t

    @reify
    def gmst(self):
        """Greenwich Mean Sidereal Time as decimal hours."""
        return sidereal_time(self)

    @reify
    def gast(self):
        """Greenwich Apparent Sidereal Time as decimal hours."""
        d_psi, _ = self._nutation_angles_radians
        tt = self.tt
        # TODO: move this into an eqeq function?
        c_terms = equation_of_the_equinoxes_complimentary_terms(tt)
        eq_eq = d_psi * cos(self._mean_obliquity_radians) + c_terms
        return (self.gmst + eq_eq / tau * 24.0) % 24.0

    # Low-precision floats generated from internal float pairs.

    @property
    def tai(self):
        return self.whole + self.tai_fraction

    @property
    def tt(self):
        return self.whole + self.tt_fraction

    @property
    def tdb(self):
        return self.whole + self.tdb_fraction

    @property
    def ut1(self):
        return self.whole + self.ut1_fraction

    # Crucial functions of time.

    def nutation_matrix(self):
        """Compute the 3×3 nutation matrix N for this date."""
        d_psi, d_eps = self._nutation_angles_radians
        mean_obliquity = self._mean_obliquity_radians
        true_obliquity = mean_obliquity + d_eps
        return build_nutation_matrix(mean_obliquity, true_obliquity, d_psi)

    def precession_matrix(self):
        """Compute the 3×3 precession matrix P for this date."""
        return compute_precession(self.tdb)

    # Various dunders.

    def __eq__(self, other_time):
        return self.__sub__(other_time) == 0.0

    def __sub__(self, other_time):
        if not isinstance(other_time, Time):
            return NotImplemented
        return ((self.whole - other_time.whole)
                + (self.tt_fraction - other_time.tt_fraction))

    def __hash__(self):
        # Someone wanted to use Time objects with functools.lru_cache so
        # we make this attempt to support hashability; beware that it
        # will return the same hash for very closely spaced times that
        # all round to the same floating point TT.
        return hash(self.tt)

    # Deprecated attributes that were once used internally, consuming
    # memory with matrices that are never used again by Skyfield once
    # t.M has been computed.

    P = reify(precession_matrix)
    N = reify(nutation_matrix)
    P.__doc__ = N.__doc__ = None  # omit from Sphinx documentation

    @reify
    def PT(self): return rollaxis(self.P, 1)
    @reify
    def NT(self): return rollaxis(self.N, 1)

def julian_day(year, month=1, day=1, julian_before=None):
    """Given a calendar date, return a Julian day integer.

    Uses the proleptic Gregorian calendar unless ``julian_before`` is
    set to a specific Julian day, in which case the Julian calendar is
    used for dates older than that.

    """
    # Support months <1 and >12 by overflowing cleanly into adjacent years.
    y, month = divmod(month - 1, 12)
    year = year + y
    month += 1

    # See the Explanatory Supplement to the Astronomical Almanac 15.11.
    janfeb = month <= 2
    g = year + 4716 - janfeb
    f = (month + 9) % 12
    e = 1461 * g // 4 + day - 1402
    J = e + (153 * f + 2) // 5

    mask = 1 if (julian_before is None) else (J >= julian_before)
    J += (38 - (g + 184) // 100 * 3 // 4) * mask
    return J

def julian_date(year, month=1, day=1, hour=0, minute=0, second=0.0):
    """Given a proleptic Gregorian calendar date and time, build a Julian date.

    The difference between a “Julian day” and a “Julian date” is that
    the “day” is the integer part, while the “date” includes a fraction
    indicating the time.

    """
    return julian_day(year, month, day) - 0.5 + (
        second + minute * 60.0 + hour * 3600.0) / DAY_S

def julian_date_of_besselian_epoch(b):
    return 2415020.31352 + (b - 1900.0) * 365.242198781

def compute_calendar_date(jd_integer, julian_before=None):
    """Convert Julian day ``jd_integer`` into a calendar (year, month, day).

    Uses the proleptic Gregorian calendar unless ``julian_before`` is
    set to a specific Julian day, in which case the Julian calendar is
    used for dates older than that.

    """
    use_gregorian = (julian_before is None) or (jd_integer >= julian_before)

    # See the Explanatory Supplement to the Astronomical Almanac 15.11.
    f = jd_integer + 1401
    f += use_gregorian * ((4 * jd_integer + 274277) // 146097 * 3 // 4 - 38)
    e = 4 * f + 3
    g = e % 1461 // 4
    h = 5 * g + 2
    day = h % 153 // 5 + 1
    month = (h // 153 + 2) % 12 + 1
    year = e // 1461 - 4716 + (12 + 2 - month) // 12
    return year, month, day

calendar_date = compute_calendar_date  # old name, in case anyone used it

def calendar_tuple(jd_float, fraction=0.0, julian_before=None):
    """Return a (year, month, day, hour, minute, second.fraction) tuple."""
    jd_float = _to_array(jd_float)
    whole1, fraction1 = divmod(jd_float, 1.0)
    whole2, fraction = divmod(fraction1 + fraction + 0.5, 1.0)
    whole = (whole1 + whole2).astype(int)
    year, month, day = compute_calendar_date(whole, julian_before)
    second = fraction * 86400.0
    minute, second = divmod(second, 60.0)
    minute = minute.astype(int)
    hour, minute = divmod(minute, 60)
    return year, month, day, hour, minute, second

def tdb_minus_tt(jd_tdb, fraction_tdb=0.0):
    """Computes how far TDB is in advance of TT, given TDB.

    Given that the two time scales never diverge by more than 2ms, TT
    can also be given as the argument to perform the conversion in the
    other direction.

    """
    t = (jd_tdb - T0 + fraction_tdb) / 36525.0

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
    ancient = load_bundled_npy('morrison_stephenson_deltat.npy')
    historic = load_bundled_npy('historic_deltat.npy')

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

_format_uses_milliseconds = re.compile(r'%[-_0^#EO]*f').search
_format_uses_seconds = re.compile(r'%[-_0^#EO]*[STXc]').search
_format_uses_minutes = re.compile(r'%[-_0^#EO]*[MR]').search
_format_uses_day_of_year = re.compile(r'%[-_0^#EO]*j').search

def _datetime_to_utc_tuple(dt):
    z = dt.tzinfo
    if z is None:
        raise ValueError(_naive_complaint)
    if z is not utc:
        dt = dt.astimezone(utc)
    return (dt.year, dt.month, dt.day,
            dt.hour, dt.minute, dt.second + dt.microsecond / 1e6)

def _normalize_jd_and_fraction(jd, fraction):
    jd = _to_array(jd)
    if fraction is None:
        jd, fraction = divmod(jd, 1.0)
    else:
        jd, fraction = _reconcile(jd, _to_array(fraction))
    return jd, fraction

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

Or use a timezone object like those provided by the third-party `pytz` library:

    from pytz import timezone
    eastern = timezone('US/Eastern')
    d = eastern.localize(datetime(2014, 1, 16, 1, 32, 9))"""
