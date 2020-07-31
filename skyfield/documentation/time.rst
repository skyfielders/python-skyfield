
================
 Dates and Time
================

.. currentmodule:: skyfield.timelib

Astronomers use several different numerical scales for measuring time.
Skyfield often has to use several timescales
even within a single computation.
So the :class:`Time` class
is designed to cache each new time scale
when a calculation first demands it.
Further demands for the same time scale
can then be satisfied without recomputing the value again.

Each time scale supported by :class:`Time`
is described in detail in one of the sections below.
The supported time scales are:

* ``t.utc`` — Coordinated Universal Time (“Greenwich Time”)
* ``t.tai`` — International Atomic Time
* ``t.tt`` — Terrestrial Time
* ``t.tdb`` — Barycentric Dynamical Time (the JPL’s *T*\ :sub:`eph`)
* ``t.ut1`` — Universal Time

To specify a time,
first build a :class:`Timescale` object
by calling Skyfield’s ``load.timescale()`` routine.
This downloads several data files from international authorities —
the United States Naval Observatory
and the International Earth Rotation Service —
to make sure that Skyfield has current information
about both leap seconds and the orientation of the Earth.
(Both topics are covered in more detail below.)

Once you have a timescale object,
which Skyfield programmers conventionally name ``ts``,
you can use its methods to create times
specified using any of the time scales listed above:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup()

   import numpy as np
   np.set_printoptions(suppress=True)

.. testcode::

    # Building a date object

    from skyfield.api import load
    ts = load.timescale()
    t = ts.utc(2014, 1, 18)

The possibilities that will be explored in the course of this page
are::

    # All the ways you can create a Time object
    # using a timescale:

    t = ts.utc(year, month, day, hour, minute, second)
    t = ts.utc(dt)        # Python datetime.datetime object

    t = ts.tai(year, month, day, hour, minute, second)
    t = ts.tai_jd(float)  # Julian date

    t = ts.tt(year, month, day, hour, minute, second)
    t = ts.tt_jd(float)   # Julian date

    t = ts.tdb(year, month, day, hour, minute, second)
    t = ts.tdb_jd(float)  # Julian date

    t = ts.ut1(year, month, day, hour, minute, second)
    t = ts.ut1_jd(float)  # Julian date

Once you have constructed a :class:`Time` object,
you can provide it to any Skyfield routine that needs it.

.. testcode::

    from skyfield.api import load

    planets = load('de421.bsp')
    earth = planets['earth']

    # Building a date and using it with at()

    ts = load.timescale()
    t = ts.utc(2014, 1, 1)
    print(earth.at(t).position.au)

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

If you will need to use the same time value several times
then it is best to create the object once,
through a single method call to your timescale object,
and then use that single time repeatedly in your calculations.
Not only will you avoid asking Skyfield to repeatedly translate
the same time value between the different time scales,
but other expensive values that depend upon time
are also automatically cached on the date object.
(See the section on :ref:`date-cache` for more details.)

.. _building-dates:

Building and printing UTC
=========================

The ``utc`` parameter in the examples above
specifies Coordinated Universal Time (UTC),
the world clock known affectionately as “Greenwich Mean Time”
which is the basis for all of the world’s timezones.

If you are comfortable dealing directly with UTC,
you can simply set and retrieve it manually.
You can provide its constructor with just the year, month, and day,
or be more specific and give an hour, minute, and second.
And not only can you attach a fraction to the seconds,
but you can also freely use fractional days, hours, and minutes.
For example:

.. testcode::

    # Four ways to specify 2014 January 18 01:35:37.5

    t1 = ts.utc(2014, 1, 18.06640625)
    t2 = ts.utc(2014, 1, 18, 1.59375)
    t3 = ts.utc(2014, 1, 18, 1, 35.625)
    t4 = ts.utc(2014, 1, 18, 1, 35, 37.5)

    assert t1 == t2 == t3 == t4    # True!

    # Several ways to print a time as UTC.

    print(tuple(t1.utc))
    print(t1.utc_iso(' '))
    print(t1.utc_jpl())
    print(t1.utc_strftime('Date %Y-%m-%d and time %H:%M:%S'))

.. testoutput::

    (2014, 1, 18, 1, 35, 37.5)
    2014-01-18 01:35:38Z
    A.D. 2014-Jan-18 01:35:37.5000 UT
    Date 2014-01-18 and time 01:35:38

The 6 values returned by ``utc()``
can also be accessed as the attributes
``year``, ``month``, ``day``, ``hour``, ``minute``,  and ``second``.

.. testcode::

    print(t1.utc.year, '/', t2.utc.month, '/', t2.utc.day)
    print(t2.utc.hour, ':', t2.utc.minute, ':', t2.utc.second)

.. testoutput::

    2014 / 1 / 18
    1 : 35 : 37.5

And by scraping together the minimal support for UTC
that exists in the Python Standard Library,
Skyfield is able to offer a :func:`~skyfield.timelib.Timescale.now()` function
that reads your system clock
and returns the current time as a Julian date object
(assuming that your operating system clock is correct
and configured with the correct time zone):

.. testcode::

    from skyfield.api import load

    # Asking the current date and time

    ts = load.timescale()
    t = ts.now()
    print(t.utc_jpl())

.. testoutput::

    A.D. 2015-Oct-11 10:00:00.0000 UT

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

But this documentation will focus on the approach
which works universally across all Python versions.
You can install the third-party `pytz`_ library
by listing it in the dependencies of your package,
or adding it to your project’s `requirements.txt`_ file,
or simply installing it manually::

    pip install pytz

Once it is installed,
building Julian dates from local times is simple.
Instantiate a normal Python ``datetime``,
pass it to the ``localize()`` method of your time zone,
and pass the result to Skyfield:

.. testcode::

    from datetime import datetime
    from pytz import timezone

    eastern = timezone('US/Eastern')

    # Converting US Eastern Time to a Julian date.

    d = datetime(2014, 1, 16, 1, 32, 9)
    e = eastern.localize(d)
    t = ts.from_datetime(e)

And if Skyfield returns a Julian date at the end of a calculation,
you can ask the Julian date object to build a ``datetime`` object
for either UTC or for your own timezone:

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
when building a Julian date.
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

    A.D. 2012-Jun-30 23:59:58.0000 UT
    A.D. 2012-Jun-30 23:59:59.0000 UT
    A.D. 2012-Jun-30 23:59:60.0000 UT
    A.D. 2012-Jul-01 00:00:00.0000 UT
    A.D. 2012-Jul-01 00:00:01.0000 UT

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
A graph will show a planet moving slightly farther
during an hour that was lengthened by a leap second.
An Earth satellite’s velocity will seem higher
when you reach the minute that includes 61 seconds.
And so forth.
Problems like these are the reason
that the :class:`Time` class only uses UTC for input and output,
and insists on keeping time internally
using the uniform time scales discussed below in :ref:`tai-tt-tdb`.

Second, leap seconds disqualify the Python ``datetime``
from use as a general way to represent time
because it refuses to accept seconds greater than 59:

.. testcode::

    datetime(2012, 6, 30, 19, 59, 60)

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: second must be in 0..59

That is why Skyfield offers a second version
of each method that returns a ``datetime``::

    t.utc_datetime_and_leap_second()
    t.astimezone_and_leap_second(tz)

These more accurate alternatives also return a ``leap_second``,
which usually has the value ``0`` but jumps to ``1``
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

.. _date-arrays:

Date arrays
===========

If you want to ask where a planet or satellite was
at a whole list of different times and dates,
then Skyfield will work most efficiently
if you build a single :class:`Time` object
that holds an entire array of dates,
instead of building many separate :class:`Time` objects.
There are three techniques for building arrays.

* Provide ``ts.utc()`` with a Python list of ``datetime`` objects.

* Provide ``tai()`` or ``tt()`` or ``tdb()`` or ``ut1()``
  with an entire NumPy array or Python list of floating point values.

* When specifying year, month, day, hour, minute, and second,
  make one of the values a list or array.

The last possibility is generally the one that is the most fun,
because its lets you vary whichever time unit you want
while holding the others steady.
And you are free to provide out-of-range values
and leave it to Skyfield to work out the correct result.
Here are some examples::

    ts.utc(range(1900, 1950))     # Fifty years 1900–1949
    ts.utc(1980, range(1, 25))    # Twenty-four months
    ts.utc(2005, 5, [1, 10, 20])  # 1st, 10th, and 20th of May

    # The ten seconds crossing the 1974 leap second
    ts.utc(1975, 1, 1, 0, 0, range(-5, 5))

The resulting :class:`Time` object will hold an array of times
instead of just a single scalar value.
As illustrated in the previous section (on leap seconds),
you can use a Python ``for`` to print each time separately:

.. testcode::

    t = ts.utc(2020, 6, 16, 7, range(4))

    for ti in t:
        print(ti.utc_strftime('%Y-%m-%d %H:%M'))

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

    t = ts.utc(2014, 1, 1)
    pos = earth.at(t).position.au
    print(pos)

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

.. testcode::

    # Whole array of Earth positions

    days = [1, 2, 3, 4]
    t = ts.utc(2014, 1, days)
    pos = earth.at(t).position.au
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

Finally, converting an array Julian date back into a calendar tuple
results in the year, month, and all of the other values
being as deep as the array itself:

.. testcode::

    print(t.utc)

.. testoutput::

    [[2014. 2014. 2014. 2014.]
     [   1.    1.    1.    1.]
     [   1.    2.    3.    4.]
     [   0.    0.    0.    0.]
     [   0.    0.    0.    0.]
     [   0.    0.    0.    0.]]

Again, simply slice across the second dimension of the array
to pull a particular calendar tuple out of the larger result:

.. testcode::

    print(t.utc[:,2])

.. testoutput::

    [2014.    1.    3.    0.    0.    0.]

The rows can be fetched not only by index
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

Did you notice how negative years work?
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

    # 2014 January 1 as a Julian Date

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
The official leap second table
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
plus exactly 32.184 seconds.
So it is now more than a minute ahead of UTC.

Barycentric Dynamical Time (TDB) runs at approximately the rate
that an atomic clock would run
if it were at rest with respect to the Solar System barycenter,
and therefore unaffected by the Earth’s motion.
The acceleration that Earth experiences in its orbit —
sometimes speeding up, sometimes slowing down —
varies the rate at which our atomic clocks
seem to run to an outside observer,
as predicted by Einstein’s theory of General Relativity.
So physical simulations of the Solar System tend to use TDB,
which is continuous with the *T*\ :sub:`eph` time scale
traditionally used for Solar System and spacecraft simulations
at the Jet Propulsion Laboratory.

UT1 and ∆T
==========

Finally, UT1 is the least uniform time scale of all
because its clock cannot be housed in a laboratory,
nor is its rate established by any human convention.
It is, rather, the clock
whose “hand” is the rotation of the Earth itself!

The UT1 time is derived from the direction
that the Earth happens to be pointing at any given moment.
And the Earth is a young world
with a still-molten iron core, a viscous mantle,
and continents that rise and fall
as each passing ice age weighs them down with ice and then melts away.
We think that we can predict, with high accuracy,
where the planets will be in their orbits
thousands of years from now.
But to predict the fluid dynamics of an elastic rotating ellipsoid
is, at the moment, beyond us.
We cannot, for example, run a simulation or formula
to predict leap seconds more than a few months ahead of time!
Instead, we simply have to watch with sensitive instruments
to see what the Earth will do next.

If you are interested in the Earth as a dynamic body,
visit the `Long-term Delta T
<http://www.usno.navy.mil/USNO/earth-orientation/eo-products/long-term>`_
page provided by the United States Naval Observatory.
You will find graphs and tables
showing how the length of Earth’s day
expands and contracts by milliseconds over the decades.
The accumulated error at any given moment
is provided as ∆T,
the evolving difference between TT and UT1
that dropped below zero in 1871 but then rose past it in 1902
and now stands at more than +67.2 seconds.

The task of governing leap seconds can be stated, then,
as the task of keeping the difference between TT and UTC
close to the natural value ∆T out in the wild.
The standards bodies promise, in fact,
that the difference between these two artificial time scales
will always be within 0.9 seconds of the observed ∆T value.

In calculations that do not involve Earth’s rotation,
∆T never arises.
The positions of planets,
the distance to the Moon,
and the movement of a comet or asteroid
all ignore ∆T completely.
When, then, does ∆T come into play?

* ∆T is used when you specify your geographic location
  as a :class:`~skyfield.toposlib.Topos`
  and Skyfield needs to compute its location at a given date and time.

* ∆T is needed to determine directions
  like “up,” “north,” and “east” when you want Skyfield
  to compute the altitude and azimuth of an object
  in your local sky.

* ∆T determines the Earth orientation for Skyfield
  when an Earth satellite position generated from TLE elements
  gets translated into a full Solar System position.

When you create your ``ts`` timescale object
at the beginning of your program,
Skyfield downloads up-to-date ``deltat.data`` and ``deltat.preds`` files
(if they are not already downloaded)
from the United States Naval Observatory.
These provide sub-millisecond level measurements
of the direction that the Earth is pointing,
allowing Skyfield to make

When you ask about dates in the far future or past,
Skyfield will run off the end of its tables
and will instead use the formula of Morrison and Stephenson (2004)
to estimate when day and night might have occurred in that era.

Setting a Custom Value For ∆T
=============================

If you ever want to specify your own value for ∆T,
then provide a ``delta_t`` keyword argument
when creating your timescale:

.. testcode::

    load.timescale(delta_t=67.2810).utc((2014, 1, 1))

.. _time-precision:

Time precision is around ~20.1 µs
=================================

Skyfield stores time vectors internally
as NumPy 64-bit floating point arrays of Julian times.
As explained in the United States Naval Observatory’s
`AA Technical Note 2011-02,
“The Error in the Double Precision Representation of Julian Dates,”
<http://aa.usno.navy.mil/software/novas/USNOAA-TN2011-02.pdf>`_
this provides fairly high precision:

    “An evaluation of the error associated with representing Julian
    dates in IEEE 754 double precision floating-point numbers
    demonstrates that Julian dates near the current epoch can be
    represented to a precision not worse than 20.1 microseconds.”

Skyfield’s own routines for turning time into strings
do careful enough rounding that you should never see effects that small.
For example, Skyfield renders the seconds of this time
attractively all the way down to 4 decimal places:

.. testcode::

    t = ts.utc(2014, 1, 18, 1, 35, 37)
    print(t.utc_jpl())

.. testoutput::

    A.D. 2014-Jan-18 01:35:37.0000 UT

It’s only if you accidentally let Python print out a raw floating point value
that you’ll see the limit of the precision:

.. testcode::

    print(t.utc.hour, t.utc.minute, t.utc.second)

.. testoutput::

    1 35 36.999999999991815

To avoid ugly output like this,
you should use Skyfield’s own time display methods
like :meth:`~skyfield.timelib.Time.utc_iso()`
and :meth:`~skyfield.timelib.Time.utc_jpl()`
or those of the ``datetime`` that Skyfield returns
from :meth:`~skyfield.timelib.Time.utc_datetime()`
when printing dates to the screen.

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
    and a vector in the dynamical reference system of this Julian date.

``C``, ``CT``
    This 3×3 matrix and its inverse
    perform the complete rotation between a vector in the ICRF
    and a vector in the celestial intermediate reference system (CIRS) of
    this Julian date.

You will typically never need to access these matrices yourself,
as they are used automatically
by the :meth:`~skyfield.positionlib.ICRF.radec()`
method when you use its  ``epoch=`` parameter
to ask for a right ascension and declination
in the dynamical reference system,
and when you ask a :class:`~skyfield.toposlib.Topos` object
for its position.

.. _matplotlib: http://matplotlib.org/
.. _pytz: http://pytz.sourceforge.net/
.. _requirements.txt: https://pip.pypa.io/en/latest/user_guide.html#requirements-files

.. testcleanup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()
