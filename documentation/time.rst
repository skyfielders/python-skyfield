
================
 Dates and Time
================

.. currentmodule:: skyfield.timelib

Astronomers use a variety of different scales to measure time.
Skyfield often has to use several timescales within a single computation!
The `Time` class is how Skyfield represents either a single moment in time
or a whole array of moments,
and keeps track of all of the different designations
assigned to that moment
by the various standard time scales:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup()

.. testcode::

    from skyfield.api import load
    ts = load.timescale()
    t = ts.tt(2000, 1, 1, 12, 0)

    print('TT date and time: ', t.tt_strftime())
    print('TAI date and time:', t.tai_strftime())
    print('UTC date and time:', t.utc_strftime())
    print('TDB Julian date: {:.10f}'.format(t.tdb))
    print('Julian century: {:.1f}'.format(t.J))

.. testoutput::

    TT date and time:  2000-01-01 12:00:00 TT
    TAI date and time: 2000-01-01 11:59:28 TAI
    UTC date and time: 2000-01-01 11:58:56 UTC
    TDB Julian date: 2451544.9999999991
    Julian century: 2000.0

The `Timescale` object returned by ``load.timescale()``
manages the conversions between different time scales
and is also how the programmer builds `Time` objects for specific dates.
Most applications create only one `Timescale` object,
which Skyfield programmers conventionally name ``ts``,
and use it to build all of their times.

For quick reference,
here are the supported timescales:

* UTC — Coordinated Universal Time (“Greenwich Time”)
* UT1 — Universal Time
* TAI — International Atomic Time
* TT — Terrestrial Time
* TDB — Barycentric Dynamical Time (the JPL’s *T*\ :sub:`eph`)

And here are links to the API documentation for time scales and times:

* :ref:`api-Timescale`
* :ref:`api-Time`

.. _choice of calendars:

Ancient and modern dates
========================

Skyfield normally uses the modern Gregorian calendar,
even for dates in history before the Gregorian calendar’s adoption in 1582.
This “proleptic” use of Gregorian dates
makes date calculations simple,
is compatible with Python’s ``datetime``,
and is also the behavior of the United States Naval Observatory library
on which many Skyfield routines were originally based.

But the Gregorian calendar is awkward
for historians and students of ancient astronomy,
because the calendar in actual use before 1582
was the old Julian calendar
established by Julius Caesar’s calendar reform in 45 BC.
The two calendars agree over the century
between the leap day of AD 200 and the leap day of AD 300.
But because the Julian calendar is not quite synchronized with the seasons,
its dates run ahead of the Gregorian calendar before that century
and run behind the Gregorian calendar after it.

If you would like Skyfield
to switch to the Julian calendar for historical dates —
both when interpreting the dates you input,
and when producing calendar dates as output —
simply give your ``Timescale`` object
the `Julian day <https://en.wikipedia.org/wiki/Julian_day>`_
on which you would like the calendar to switch.

.. testcode::

    from skyfield.api import GREGORIAN_START

    ts.julian_calendar_cutoff = GREGORIAN_START

    t = ts.tt_jd(range(2299159, 2299163))
    for s in t.tt_strftime():
        print(s)

.. testoutput::

    1582-10-03 12:00:00 TT
    1582-10-04 12:00:00 TT
    1582-10-15 12:00:00 TT
    1582-10-16 12:00:00 TT

As you can see from these four successive days in history,
Pope Gregory had the calendar jump
directly from the old Julian calendar date 1582 October 4
to the new Gregorian calendar date 1582 October 15,
bringing the date of Easter back into sync with the equinox.
Skyfield provides two constants for popular cutoff dates:

* ``GREGORIAN_START`` — Julian day 2299161,
  on which the new Gregorian calendar went into effect in Rome.

* ``GREGORIAN_START_ENGLAND`` — Julian day 2361222,
  on which the new Gregorian calendar went into effect in England in 1752
  (the reform having initially been rejected by the English bishops,
  “Seeing that the Bishop of Rome is Antichrist,
  therefore we may not communicate with him in any thing”).

You are free to choose your own cutoff Julian day
if you are studying astronomy records from a country
that adopted the Gregorian calendar on some other date.
Russia, for example, did not adopt it until the twentieth century.
The default value,
that asks the timescale to always use Gregorian dates,
is ``None``:

.. testcode::

    ts.julian_calendar_cutoff = None

Note that even the Julian calendar becomes anachronistic
before its adoption in 45 BC,
so all dates generated by Skyfield are “proleptic” before that date.
And, of course, the Julian calendar
was local to the civilization that ringed the Mediterranean.
If you are interested in relating astronomical events
to more ancient Roman calendars,
or the calendars of other civilizations,
try searching for a third-party Python package
that supports the calendar you are interested in.

.. _building-dates:

Building and printing UTC
=========================

The ``utc`` parameter in the examples above
specifies Coordinated Universal Time (UTC),
the world clock known affectionately as “Greenwich Mean Time”
which is the basis for all of the world’s timezones.
If you are comfortable dealing directly with UTC
instead of your local timezone,
you can build and display dates
without needing any other library besides Skyfield.

You can build a `Time` from a calendar date and UTC time
using `Timescale.utc`.
Provide only as many parameters as you want —
year, month, day, hour, minute, and second —
and Skyfield will fill in the rest
by defaulting to January first and zero hours, minutes, and seconds.

Feel free to use fractional days, hours, and minutes.
Here are several ways to specify the exact same time and date:

.. testcode::

    # Four ways to specify 2014 January 18 01:35:37.5

    t1 = ts.utc(2014, 1, 18.06640625)
    t2 = ts.utc(2014, 1, 18, 1.59375)
    t3 = ts.utc(2014, 1, 18, 1, 35.625)
    t4 = ts.utc(2014, 1, 18, 1, 35, 37.5)

    assert t1 == t2 == t3 == t4    # True!

    # Several ways to print a time as UTC.

    print(tuple(t1.utc))
    print(t1.utc_iso())
    print(t1.utc_strftime())
    print(t1.utc_strftime('On %Y %b %d at %H:%M:%S'))
    print(t1.utc_jpl())

.. testoutput::

    (2014, 1, 18, 1, 35, 37.5)
    2014-01-18T01:35:38Z
    2014-01-18 01:35:38 UTC
    On 2014 Jan 18 at 01:35:38
    A.D. 2014-Jan-18 01:35:37.5000 UTC

The 6 values returned by ``utc()``
can be accessed as the attributes
``year``, ``month``, ``day``, ``hour``, ``minute``,  and ``second``.

.. testcode::

    print(t1.utc.year, '/', t1.utc.month, '/', t1.utc.day)
    print(t1.utc.hour, ':', t1.utc.minute, ':', t1.utc.second)

.. testoutput::

    2014 / 1 / 18
    1 : 35 : 37.5

If you want to use the current time,
Skyfield leverages the minimal support for UTC
in the Python Standard Library
to offer a :func:`~skyfield.timelib.Timescale.now()` function
that reads your system clock
and returns the current time as a `Time` object
(assuming that your operating system clock is correct
and configured with the correct time zone):

.. testcode::

    from skyfield.api import load

    # Asking the current date and time

    ts = load.timescale()
    t = ts.now()
    print(t.utc_jpl())

.. testoutput::

    A.D. 2015-Oct-11 10:00:00.0000 UTC

UTC and your timezone
=====================

To move beyond UTC and work with other world timezones,
you will need to install a time zone database
for your version of Python.

* Every version of Python that Skyfield supports
  will work with the `pytz`_ package described in this section.

* Python 3.6 upgraded the Standard Library ``datetime`` type
  so that the contortions of `pytz`_ are no longer necessary,
  and instead recommends
  `dateutil <https://dateutil.readthedocs.io/en/stable/>`_
  for working with timezones.
  Consult its documentation if you are interested in using it.

* Python 3.9 will offer a native
  `zoneinfo <https://docs.python.org/3.9/library/zoneinfo.html>`_
  module that for the first time brings timezone support
  into the Python Standard Library.

This documentation will focus on the first approach,
which works universally across all Python versions.
You can install the third-party `pytz`_ library
by listing it in the dependencies of your package,
adding it to your project’s `requirements.txt`_ file,
or simply installing it manually::

    pip install pytz

Once it is installed,
building a `Time` from a local time is simple.
Instantiate a normal Python ``datetime``,
pass it to the ``localize()`` method of your time zone,
and pass the result to Skyfield:

.. testcode::

    from datetime import datetime
    from pytz import timezone

    eastern = timezone('US/Eastern')

    # Converting US Eastern Time to a Skyfield Time.

    d = datetime(2014, 1, 16, 1, 32, 9)
    e = eastern.localize(d)
    t = ts.from_datetime(e)

When Skyfield returns a `Time` at the end of a calculation,
you can ask for either a UTC ``datetime``
or a ``datetime`` in your own timezone:

.. testcode::

    # UTC datetime

    dt = t.utc_datetime()
    print('UTC: ' + str(dt))

    # Converting back to an Eastern Time datetime.

    dt = t.astimezone(eastern)
    print('EST: ' + str(dt))

.. testoutput::

    UTC: 2014-01-16 06:32:09+00:00
    EST: 2014-01-16 01:32:09-05:00

As we would expect,
1:32 AM in the Eastern time zone in January
is 6:32 AM local time in Greenwich, England,
five hours to the east across the Atlantic.

Note that Skyfield’s :meth:`~Time.astimezone()` method
will detect that you are using a ``pytz`` timezone
and automatically call its ``normalize()`` method for you —
which makes sure that daylight savings time is handled correctly —
to spare you from having to make the call yourself.

If you want a :class:`Time` to hold an entire array of dates,
as discussed below in :ref:`date-arrays`,
then you can provide a list of ``datetime`` objects
to the `Timescale.from_datetimes()` method.
The UTC methods will then return whole lists of values.

.. _leap-seconds:

UTC and leap seconds
====================

The rate of Earth’s rotation is gradually slowing down.
Since the UTC standard specifies a fixed length for the second,
promises a day of 24 hours, and limits an hour to 60 minutes,
the only way to stay within the rules
while keeping UTC synchronized with the Earth
is to occasionally add an extra leap second
to one of the year’s minutes.

See :ref:`the-leap-second-table` if you are interested
in printing Skyfield’s full list of leap seconds.

The `International Earth Rotation Service <http://hpiers.obspm.fr/>`_
currently restricts itself to appending a leap second
to the last minute of June or the last minute of December.
When a leap second is inserted,
its minute counts 61 seconds numbered 00–60
instead of staying within the usual range 00–59.
One recent leap second was in June 2012:

.. testcode::

    # Display 5 seconds around a leap second

    five_seconds = [58, 59, 60, 61, 62]
    t = ts.utc(2012, 6, 30, 23, 59, five_seconds)

    for string in t.utc_jpl():
        print(string)

.. testoutput::

    A.D. 2012-Jun-30 23:59:58.0000 UTC
    A.D. 2012-Jun-30 23:59:59.0000 UTC
    A.D. 2012-Jun-30 23:59:60.0000 UTC
    A.D. 2012-Jul-01 00:00:00.0000 UTC
    A.D. 2012-Jul-01 00:00:01.0000 UTC

Note that Skyfield has no problem with a calendar tuple
that has hours, minutes, or — as in this case —
seconds that are out of range.
When we provided a range of numbers 58 through 62 as seconds,
Skyfield added exactly the number of seconds we specified
to the end of June
and let the value overflow cleanly into the beginning of July.

Keep two consequences in mind when using UTC in your calculations.

First, expect an occasional jump or discrepancy
if you are striding forward through time
using the UTC minute, hour, or day.
For example,
an hourly plot of planet’s position
will show the planet moving slightly farther
during an hour that was lengthened by a leap second
than during other hours of the year.
An Earth satellite’s velocity will seem higher
when you reach the minute that includes 61 seconds.
And so forth.
Problems like these are the reason
that the :class:`Time` class only uses UTC for input and output,
and insists on keeping time internally
using the uniform time scales discussed below in :ref:`tai-tt-tdb`.

Second, leap seconds disqualify the Python ``datetime``
from use as a general way to represent time
because in many versions of Python
the ``datetime`` refuses to accept seconds greater than 59:

.. testcode::

    datetime(2012, 6, 30, 19, 59, 60)

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: second must be in 0..59

That limitation is why Skyfield offers a second version
of each method that returns a ``datetime``.
These fancier methods return a leap-second flag as an additional return value:

* `Time.utc_datetime_and_leap_second()`
* `Time.astimezone_and_leap_second()`

The leap-second return value is usually ``0`` but jumps to ``1``
when Skyfield is forced to represent a leap second
as a ``datetime`` with the incorrect time 23:59:59.

.. testcode::

    # Asking for the leap_second flag to learn the whole story

    dt, leap_second = t.astimezone_and_leap_second(eastern)

    for dt_i, leap_second_i in zip(dt, leap_second):
        print('{0}  leap_second = {1}'.format(dt_i, leap_second_i))

.. testoutput::

    2012-06-30 19:59:58-04:00  leap_second = 0
    2012-06-30 19:59:59-04:00  leap_second = 0
    2012-06-30 19:59:59-04:00  leap_second = 1
    2012-06-30 20:00:00-04:00  leap_second = 0
    2012-06-30 20:00:01-04:00  leap_second = 0

Using calendar tuples to represent UTC times is more elegant
than using Python ``datetime`` objects
because leap seconds can be represented accurately.
If your application cannot avoid using ``datetime`` objects,
then you will have to decide
whether to simply ignore the ``leap_second`` value
or to somehow output the leap second information.

Date arithmetic
===============

Dates support a few simple math operations:

.. testcode::

    from datetime import timedelta

    t - 10                   # 10 days earlier
    t + 0.25                 # 6 hours later
    t + timedelta(hours=12)  # 12 hours later

    t2 - t1  # difference between times, in days

Raw numbers,
like ``10`` and ``0.25`` above,
specify days of Terrestrial Time —
units of exactly 24 hours of 60 minutes of 60 SI seconds,
measured in the Earth’s relativistic frame of reference.
If you increment or decrement a date
across a :ref:`leap second <leap-seconds>`,
you will notice
that the clock time returned by Skyfield’s UTC functions
is one second earlier or later than you expect.

.. _date-arrays:

Date arrays
===========

If you want to ask where a planet or satellite was
across a whole series of times and dates,
then Skyfield will work most efficiently if,
instead of building many separate :class:`Time` objects,
you build a single :class:`Time` object that holds the entire array of dates.

There are three techniques for building a `Time` array.

* Provide :meth:`~Timescale.tai()` or :meth:`~Timescale.tt()`
  or :meth:`~Timescale.tdb()` or :meth:`~Timescale.ut1()`
  with a Python list or NumPy array of numbers
  for one of the six components of the calendar date
  (year, month, day, hour, minute, or second).

* Provide :meth:`~Timescale.tai_jd()` or :meth:`~Timescale.tt_jd()`
  or :meth:`~Timescale.tdb_jd()` or :meth:`~Timescale.ut1_jd()`
  with a list or NumPy array of floating point numbers.

* Provide :meth:`~Timescale.from_datetimes()`
  with a Python list of ``datetime`` objects.

The first possibility is generally the one that is the most fun,
because its lets you vary whichever time unit you want
while holding the others constant.
You are free to provide out-of-range values
and leave it to Skyfield to work out the correct result.
Here are some examples::

    ts.utc(range(1900, 1950))     # Fifty years 1900–1949
    ts.utc(1980, range(1, 25))    # 24 months of 1980 and 1981
    ts.utc(2005, 5, [1, 11, 21])  # 1st, 11th, and 21st of May

    # Negative values work too!  Here are the
    # ten seconds crossing the 1974 leap second.
    ts.utc(1975, 1, 1, 0, 0, range(-5, 5))

The resulting :class:`Time` object will hold an array of times.
As illustrated in the previous section (on leap seconds),
you can use a Python ``for`` to print each time separately:

.. testcode::

    t = ts.utc(2020, 6, 16, 7, range(4))

    for s in t.utc_strftime('%Y-%m-%d %H:%M'):
        print(s)

.. testoutput::

    2020-06-16 07:00
    2020-06-16 07:01
    2020-06-16 07:02
    2020-06-16 07:03

When you provide a time array as input to a Skyfield calculation,
the output array will have an extra dimension
that expands what would normally be a single result
into as many results as you provided dates.
We can compute the position of the Earth as an example:

.. testcode::

    # Single Earth position

    planets = load('de421.bsp')
    earth = planets['earth']

    t = ts.utc(2014, 1, 1)
    pos = earth.at(t).xyz.au
    print(pos)

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

.. testcode::

    # Whole array of Earth positions

    days = [1, 2, 3, 4]
    t = ts.utc(2014, 1, days)
    pos = earth.at(t).xyz.au
    print(pos)

.. testoutput::

    [[-0.17461758 -0.19179872 -0.20891924 -0.22597338]
     [ 0.88567056  0.88265548  0.87936337  0.87579547]
     [ 0.38384886  0.38254134  0.38111391  0.37956709]]

Note the shape of the resulting NumPy array.
If you unpack this array into three names,
then you get three four-element arrays
corresponding to the four dates.
These four-element arrays are ready to be submitted to `matplotlib`_
and other scientific Python tools:

.. testsetup::

    def plot(*args): pass

.. testcode::

    x, y, z = pos    # four values each
    plot(x, y)       # example matplotlib call

If you instead slice along the second axis,
then you can retrieve an individual position for a particular date —
and the first position is exactly what was returned above
when we computed the January 1st position by itself:

.. testcode::

    print(pos[:,0])

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

You can combine a Python ``for`` loop with Python’s ``zip()`` builtin
to print each time alongside the corresponding coordinates.
There are two techniques,
one of which is less efficient and the other more efficient.

.. testcode::

    # Less efficient: loop over `t`, forcing the creation of
    # a separate `Time` object for each iteration of the loop.

    for ti, xi, yi, zi in zip(t, x, y, z):
        print('{}  x = {:.2f} y = {:.2f} z = {:.2f}'.format(
            ti.utc_strftime('%Y-%m-%d'), xi, yi, zi,
        ))

.. testoutput::

    2014-01-01  x = -0.17 y = 0.89 z = 0.38
    2014-01-02  x = -0.19 y = 0.88 z = 0.38
    2014-01-03  x = -0.21 y = 0.88 z = 0.38
    2014-01-04  x = -0.23 y = 0.88 z = 0.38

.. testcode::

    # More efficient: loop over the output of a `Time` method,
    # which returns an array of the same length as `t`.

    t_strings = t.utc_strftime('%Y-%m-%d')

    for tstr, xi, yi, zi in zip(t_strings, x, y, z):
        print('{}  x = {:.2f} y = {:.2f} z = {:.2f}'.format(
            tstr, xi, yi, zi,
        ))

.. testoutput::

    2014-01-01  x = -0.17 y = 0.89 z = 0.38
    2014-01-02  x = -0.19 y = 0.88 z = 0.38
    2014-01-03  x = -0.21 y = 0.88 z = 0.38
    2014-01-04  x = -0.23 y = 0.88 z = 0.38

Finally, converting an array `Time` back into a calendar tuple
results in the year, month, day, hour, minute, and second
each having the same dimension as the array itself:

.. testcode::

    print(t.utc)

.. testoutput::

    [[2014. 2014. 2014. 2014.]
     [   1.    1.    1.    1.]
     [   1.    2.    3.    4.]
     [   0.    0.    0.    0.]
     [   0.    0.    0.    0.]
     [   0.    0.    0.    0.]]

Simply slice across the second dimension of the array
to pull a particular calendar tuple out of the larger result:

.. testcode::

    print(t.utc[:,2])

.. testoutput::

    [2014.    1.    3.    0.    0.    0.]

Slicing in the other direction,
the rows can be fetched not only by index
but also through the attribute names ``year``, ``month``, ``day``,
``hour``, ``minute``,  and ``second``.

.. testcode::

    print(t.utc.year)
    print(t.utc.month)
    print(t.utc.day)
    print(t.utc.hour)

.. testoutput::

    [2014. 2014. 2014. 2014.]
    [1. 1. 1. 1.]
    [1. 2. 3. 4.]
    [0. 0. 0. 0.]

.. _tai-tt-tdb:

Uniform time scales: TAI, TT, and TDB
=====================================

Date arithmetic becomes very simple
as we leave UTC behind and consider completely uniform time scales.
Days are always 24 hours, hours always 60 minutes,
and minutes always 60 seconds without any variation or exceptions.
Such time scales are not appropriate for your morning alarm clock
because they are never delayed or adjusted
to stay in sync with the slowing rotation of the earth.
But that is what makes them useful for astronomical calculation —
because physics keeps up its dance,
and the stars and planets move in their courses,
whether humanity pauses to observe a UTC leap second or not.

Because they make every day the same length,
uniform time scales can express dates
as a simple floating-point count of days elapsed.
To make all historical dates come out as positive numbers,
astronomers traditionally assign each date a “Julian day” number
that starts counting at 4713 BC January 1 in the old Julian calendar —
the same date as 4714 BC November 24 in our Gregorian calendar.
Following a tradition going back to the Greeks and Ptolemy,
the count starts at noon,
since the sun’s transit is an observable event
but the moment of midnight is not.

So twelve noon was the moment of Julian date zero:

.. testcode::

    # When was Julian date zero?

    bc_4714 = -4713
    t = ts.tt(bc_4714, 11, 24, 12)
    print(t.tt)

.. testoutput::

    0.0

Did you notice how negative years work —
that we expressed 4714 BC using the negative number ``-4713``?
People still counted by starting at one, not zero,
when the scholar Dionysius Exiguus created the eras BC and AD
in around the year AD 500.
So his scheme has 1 BC followed immediately by AD 1 without a break.
To avoid an off-by-one error,
astronomers usually ignore BC and count backwards through a year zero
and on into negative years.
So negative year *−n* is what might otherwise be called
either “*n+1* BC” or “*n+1* BCE” in a history textbook.

More than two million days have passed since 4714 BC,
so modern dates tend to be rather large numbers:

.. testcode::

    # 2014 January 1 00:00 UTC expressed as Julian dates

    t = ts.utc(2014, 1, 1)
    print('TAI = %r' % t.tai)
    print('TT  = %r' % t.tt)
    print('TDB = %r' % t.tdb)

.. testoutput::

    TAI = 2456658.5004050927
    TT  = 2456658.5007775924
    TDB = 2456658.500777592

What are these three different uniform time scales?

International Atomic Time (TAI) is maintained
by the worldwide network of atomic clocks
referenced by researchers with a need for very accurate time.
The official :ref:`leap second table <the-leap-second-table>`
is actually a table of offsets between TAI and UTC.
At the end of June 2012, for example,
the TAI−UTC offset was changed from 34.0 to 35.0
which is what generated the leap second in UTC.

Terrestrial Time (TT) differs from TAI
only because astronomers
were already maintaining a uniform time scale of their own
before TAI was established,
using a slightly different starting point for the day.
For practical purposes, TT is simply TAI
plus exactly 32.184 seconds.
So it is now more than a minute ahead of UTC.

You can not only ask Skyfield for TT as a Julian date and a calendar date,
but as a floating-point number of years
of exactly 365.25 days each —
a value which is often used as the time parameter
in long-term astronomical formulae:

.. testcode::

    print('Julian year = {:.4f}'.format(t.J))

.. testoutput::

    Julian year = 2014.0000

Finally,
Barycentric Dynamical Time (TDB) runs at approximately the rate
that an atomic clock would run
if it were at rest with respect to the Solar System barycenter,
and therefore unaffected by the Earth’s motion.
The acceleration that Earth experiences in its orbit —
sometimes speeding up, sometimes slowing down —
varies the rate at which our atomic clocks
run relative to an outside observer,
as predicted by Einstein’s theory of General Relativity.
So physical simulations of the Solar System use TDB as their clock.
It is considered equivalent to the *T*\ :sub:`eph` time scale
traditionally used for Solar System and spacecraft simulations
at the Jet Propulsion Laboratory.

.. _downloading-timescale-files:

UT1 and downloading IERS data
=============================

Finally, UT1 is the least uniform time scale of all
because its clock cannot be housed in a laboratory,
nor is its rate established by any human convention.
It is, rather, the clock
whose “hand” is the rotation of the Earth itself!
The direction that the Earth is facing determines
not only the coordinates of every city and observatory in the world,
but also the local directions that each site
will designate as their local “up”, “north”, and “east”.

It is hard to predict future values for UT1.
The Earth is a young world
with a still-molten iron core,
a viscous mantle,
and ice ages that move water weight into glaciers at the poles
then release it back into the ocean.
While we think we can predict, for example,
Jupiter’s position thousands of years from now,
predicting the fluid dynamics of the elastic rotating ellipsoid we call home
is — at the moment — beyond us.
We can only watch with sensitive instruments
to see what the Earth does next.

Skyfield relies on the IERS,
the International Earth Rotation Service,
for accurate measurements of UT1
and for the schedule of leap seconds (discussed above)
that keeps UTC from straying more than 0.9 seconds away from UT1.

Each new version of Skyfield carries recent IERS data in internal tables.
This data will gradually fall out of date after each Skyfield release,
however, with two consequences:

* The next time the IERS declares a new leap second
  that is not listed in Skyfield’s built-in tables,
  Skyfield’s UTC time will be off by 1 second
  for every date that falls after the leap second.

* As the Earth’s rotation speeds up or slows down in the coming years
  more than was predicted in Skyfield’s built-in UT1 tables,
  Skyfield’s idea of where the Earth is pointing will grow less accurate.
  This will affect both the position and direction
  of each :class:`~skyfield.toposlib.GeographicPosition` —
  whether used as an observer or a target —
  and will also affect Earth satellite positions.

You can avoid both of these problems
by periodically downloading new data from the IERS.
Simply specify that you don’t want Skyfield to use its builtin tables.
In that case :meth:`~skyfield.iokit.Loader.timescale()`
will instead download ``finals2000A.all`` from the IERS:

::

    # Download and use the `finals.all` file.

    ts = load.timescale(builtin=False)

::

    [#################################] 100% finals2000A.all

As usual with data files,
Skyfield will only download the file the first time you need it,
then will keep using that same copy of the file that it finds on disk.

Note that the international agencies responsible for the file’s distribution
sometimes have trouble keeping their servers up.
For example, as I write this in May of 2022,
the file cannot be fetched from ``ftp.iers.org``
because of an
`Outage of iers.org data servers <https://www.iers.org/IERS/EN/NewsMeetings/News/news_001.html>`_
reported on their website.
At
`Skyfield issue #730 <https://github.com/skyfielders/python-skyfield/issues/730>`_
and
`Skyfield issue #732 <https://github.com/skyfielders/python-skyfield/issues/732>`_
you can find links to alternative data sources
which various Skyfield users have been able to access in the meantime.

If your script will always have Internet access
and you worry about the file falling out of date
(and if you can trust the “modify time” file attribute on your filesystem),
then you can have Skyfield download a new copy
once the file on disk has grown too old
(where “too old” for your application
must be determined by comparing your accuracy needs
with how quickly UT1 diverges without fresh IERS data;
this example uses 30 days only as an illustration):

::

    if load.days_old('finals2000A.all') > 30.0:
        load.download('finals2000A.all')

    ts = load.timescale(builtin=False)

But, beware!
For compatibility with versions of Skyfield ≤ 1.30,
Skyfield will ignore ``finals2000A.all``
if the three old files
``deltat.data``, ``deltat.preds``, and ``Leap_Second.dat``
exist in the loader’s directory,
in which case it will use them instead.
This is to prevent users who specify ``builtins=False``,
but who downloaded the three necessary files long ago,
from experiencing an unexpected download attempt.
The hope is that all scripts
which did not previously need Internet access
will continue to run without it.

If you ever want to display or plot the behavior of UT1,
astronomers use two common conventions
for stating the difference between clock time and UT1.
Skyfield supports them both.

.. testcode::

    print('{:+.4f}'.format(t.dut1))
    print('{:+.4f}'.format(t.delta_t))

.. testoutput::

    -0.0970
    +67.2810

The two quantities are:

* DUT1 — The difference between UTC and UT1,
  which should always be less than 0.9 if the IERS succeeds at its mission.
  Note that there are two different reasons that this value changes:
  every day it changes a small amount because of the drift of UT1;
  and superimposed on this drift
  is a big jump of 1.0 seconds
  every time a leap second passes.

* ∆T — The difference between TT and UT1.
  This is a much more straightforward value than DUT1,
  without all the ugly discontinuities caused by leap seconds.
  Because TT is a uniform timescale,
  ∆T provides a pure and continuous record
  of how UT1 has changed over the decades
  that we have been measuring the Earth’s rotation to high precision.

.. _custom-delta-t:

Setting a Custom Value For ∆T
=============================

If you ever want to specify your own value for ∆T,
then provide a ``delta_t`` keyword argument
when creating your timescale:

.. testcode::

    load.timescale(delta_t=67.2810).utc((2014, 1, 1))

.. _date-cache:

Values cached on the Time object
================================

When you create a :class:`Time`
it goes ahead and computes its ``tt`` Terrestrial Time attribute
starting from whatever time argument you provide.
If you provide the ``utc`` parameter, for example,
then the date first computes and sets ``tai``
and then computes and sets ``tt``.
Each of the other time attributes only gets computed once,
the first time you access it.

The general rule is that attributes are only computed once,
and can be accessed again and again for free,
while methods never cache their results —
think of the ``()`` parentheses after a method name
as your reminder that “this will do a fresh computation every time.”

In addition to time scales,
each :class:`Time` object caches several other quantities
that are often needed in astronomy.
Skyfield only computes these attributes on-demand,
the first time the user tries to access them
or invokes a computation that needs their value:

``gmst``
    Greenwich Mean Sidereal Time in hours,
    in the range 0.0 ≤ ``gmst`` < 24.0.

``gast``
    Greenwich Apparent Sidereal Time in hours,
    in the range 0.0 ≤ ``gast`` < 24.0.

``M``, ``MT``
    This 3×3 matrix and its inverse
    perform the complete rotation between a vector in the ICRF
    and a vector in the dynamical reference system for this time and date.

``C``, ``CT``
    This 3×3 matrix and its inverse
    perform the complete rotation between a vector in the ICRF
    and a vector in the celestial intermediate reference system (CIRS)
    of this time and date.

You will typically never need to access these matrices yourself,
as they are used automatically
by the :meth:`~skyfield.positionlib.ICRF.radec()`
method when you use its  ``epoch=`` parameter
to ask for a right ascension and declination
in the dynamical reference system,
and when you ask a :class:`~skyfield.toposlib.GeographicPosition` object
for its position.

.. _matplotlib: http://matplotlib.org/
.. _pytz: http://pytz.sourceforge.net/
.. _requirements.txt: https://pip.pypa.io/en/latest/user_guide.html#requirements-files

.. testcleanup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()
