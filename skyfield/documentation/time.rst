
================
 Dates and Time
================

.. currentmodule:: skyfield.api

Astronomers use several different numerical scales for measuring time.
Since Skyfield often has to reference several such scales
over the course of a single calculation,
the :class:`JulianDate` class
is designed to cache each new time scale
when a calculation first demands it.
Further demands for the same time scale
can then be satisfied without recomputing the value again.

Each time scale supported by :class:`JulianDate`
is described in detail in one of the sections below.
The supported time scales are:

* ``jd.utc`` — Coordinated Universal Time (“Greenwich Time”)
* ``jd.tai`` — International Atomic Time
* ``jd.tt`` — Terrestrial Time
* ``jd.tdb`` — Barycentric Dynamical Time (the JPL’s “T\ :sub:`eph`”)
* ``jd.ut1`` — Universal Time

To build a Julian date object,
simply use one of these time scales
as a keyword argument to specify the moment you want to represent.

.. testcode::

    from skyfield.api import JulianDate
    JulianDate(utc=(2014, 1, 18))

The possibilities that will be explored in the course of this page
are::

    Parameter to the JulianDate constructor
    │
    ├── utc = Calendar tuple or Python `datetime` or `date`
    ├── tai = Calendar tuple or Julian date float
    ├── tt  = Calendar tuple or Julian date float
    ├── tdb = Calendar tuple or Julian date float
    └── ut1 = Calendar tuple or Julian date float

    Calendar tuple
    │
    ├── (year, month, day)
    ├── (year, month, day, hour)
    ├── (year, month, day, hour, minute)
    └── (year, month, day, hour, minute, second)

There are two ways to provide a date argument
to a Skyfield routine that needs one.
You can either supply a :class:`JulianDate`
that you have gone ahead and built yourself,
or you can simply call the routine with the same keyword argument
that you would have used to build the Julian date
and Skyfield will construct the date internally.
So the following two calculations for January 1st
come out exactly the same:

.. testcode::

    from skyfield.api import load
    planets = load('de421.bsp')
    earth = planets['earth']

    jd = JulianDate(utc=(2014, 1, 1))
    print(earth.at(jd).position.au)

    print(earth.at(utc=(2014, 1, 1)).position.au)

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]
    [-0.17461758  0.88567056  0.38384886]

Be sure to build the Julian date yourself
if you have the opportunity to use it over and over again
in several different calculations.
Not only will you avoid recomputing the time scales,
but other expensive calculations will be cached on the date object
(see the section on :ref:`date-cache` for more details).

.. _building-dates:

UTC and your timezone
=====================

The ``utc`` parameter in the examples above
specifies Coordinated Universal Time (UTC),
the world clock known affectionately as “Greenwich Mean Time”
which is the basis for all of the world’s timezones.

If you are comfortable dealing directly with UTC,
you can simply set and retrieve it manually.
A Python tuple is the most convenient way
to represent the year, month, and day of a calendar date:

.. testcode::

    # Four ways to specify 2014 January 18 01:35:37.5

    jd  = JulianDate(utc=(2014, 1, 18.06640625))
    jd2 = JulianDate(utc=(2014, 1, 18, 1.59375))
    jd3 = JulianDate(utc=(2014, 1, 18, 1, 35.625))
    jd4 = JulianDate(utc=(2014, 1, 18, 1, 35, 37.5))

    assert jd == jd2 == jd3 == jd4

    # Retrieving UTC as either a tuple or string

    print(jd.utc)
    print(jd.utc_iso())
    print(jd.utc_jpl())
    print(jd.utc_strftime('Date %Y-%m-%d and time %H:%M:%S'))

.. testoutput::

    (2014, 1, 18, 1, 35, 37.5)
    2014-01-18T01:35:38Z
    A.D. 2014-Jan-18 01:35:37.5000 UT
    Date 2014-01-18 and time 01:35:38

And by scraping together the minimal support for UTC
that exists in the Python Standard Library,
Skyfield is able to offer a :func:`~skyfield.api.now()` function
that reads your system clock
and returns the current time as a Julian date object
(assuming that your operating system clock is correct
and configured with the correct time zone):

.. testsetup::

    import numpy as np
    np.set_printoptions(suppress=True)

    from skyfield import api
    def now():
        """Return a constant "now"."""
        return api.JulianDate(utc=(2014, 1, 18, 23, 10, 9))
    api.now = now

.. testcode::

    from skyfield.api import now

    jd = now()
    print(jd.utc_jpl())

.. testoutput::

    A.D. 2014-Jan-18 23:10:09.0000 UT

But to move beyond UTC to working with actual timezones,
you will need to install
the third-party `pytz`_ package,
either by listing it in the dependencies of your package,
adding it to your project’s `requirements.txt`_ file,
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
    jd = JulianDate(utc=e)

And if Skyfield returns a Julian date at the end of a calculation,
you can ask the Julian date object to build a ``datetime`` object
for either the UTC date or your own timezone:

.. testcode::

    # UTC datetime

    dt = jd.utc_datetime()
    print('UTC: ' + str(dt))

    # Converting back to an Eastern Time datetime.

    dt = jd.astimezone(eastern)
    print('EST: ' + str(dt))

.. testoutput::

    UTC: 2014-01-16 06:32:09+00:00
    EST: 2014-01-16 01:32:09-05:00

As we would expect,
1:32 AM in the Eastern time zone in January
is 6:32 AM local time in Greenwich, England,
five hours to the east across the Atlantic.

Note that Skyfield’s :meth:`~JulianDate.astimezone()` method
will detect that you are using a ``pytz`` timezone
and automatically call its ``normalize()`` method for you —
which makes sure that daylight savings time is handled correctly —
to spare you from having to make the call yourself.

If you want a :class:`JulianDate` to hold an entire array of dates,
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
The most recent leap second was in June 2012:

.. testcode::

    five_seconds = range(58, 58 + 5)
    tup = (2012, 6, 30, 23, 59, five_seconds)
    jd = JulianDate(utc=tup)

    for string in jd.utc_jpl():
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
It simply adds as many seconds as we specified
and lets the value overflow cleanly into the beginning of July.

Keep two consequences in mind when using UTC in your calculations.

First, expect an occasional jump or discrepancy
if you are striding forward through time
using the UTC minute, hour, or day.
A graph will show a planet moving slightly farther
during an hour that was lengthened by a leap second;
an Earth satellite’s velocity will seem higher
when you reach the minute that includes 61 seconds;
and so forth.
Problems like these are the reason
that the :class:`JulianDate` only uses UTC for input and output,
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

    jd.utc_datetime_and_leap_second()
    jd.astimezone_and_leap_second(tz)

These more accurate alternatives also return a ``leap_second``,
which usually has the value ``0`` but jumps to ``1``
when Skyfield is forced to represent a leap second
as a ``datetime`` with the incorrect time 23:59:59.

.. testcode::

    dt, leap_second = jd.astimezone_and_leap_second(eastern)

    for dt_i, leap_second_i in zip(dt, leap_second):
        print(str(dt_i) + ' leap_second = ' + str(leap_second_i))

.. testoutput::

    2012-06-30 19:59:58-04:00 leap_second = 0
    2012-06-30 19:59:59-04:00 leap_second = 0
    2012-06-30 19:59:59-04:00 leap_second = 1
    2012-06-30 20:00:00-04:00 leap_second = 0
    2012-06-30 20:00:01-04:00 leap_second = 0

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

Skyfield works most efficiently
when you build a single :class:`JulianDate` object
that holds an entire array of dates,
instead of building many separate :class:`JulianDate` objects.
There are three techniques for building arrays.

* Make ``utc=`` a list of ``datetime`` objects.

* Specify ``tai=`` or ``tt=`` or ``tdb=`` or ``ut1=``
  using an entire NumPy array or Python list of floating point values.

* With any parameter,
  use a calendar tuple with one element
  set to a whole list or array of values
  instead of just being a single value.

The last possibility is generally the one that is the most fun,
because its lets you vary whichever time unit you want
while holding the others steady.
And you are free to provide out-of-range values
and leave it to Skyfield to work out the correct result.
Here are some examples::

    utc=(range(1900, 1950),)    # Fifty years
    utc=(1980, range(1, 25))    # Twenty-four months
    utc=(2005, 5, [1, 10, 20])  # 1st, 10th, 20th of May

    # The ten seconds crossing the 1974 leap second
    utc=(1975, 1, 1, 0, 0, range(-5, 5))

When you provide an array :class:`JulianDate` to a Skyfield calculation,
the resulting array will have an extra dimension
expanding what would normally be a single result
into as many results as you provided dates.
We can compute the position of the Earth as an example:

.. testcode::

    # Single Earth position

    jd = JulianDate(utc=(2014, 1, 1))
    pos = earth.at(jd).position.au
    print(pos)

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

.. testcode::

    # Whole array of Earth positions

    days = [1, 2, 3, 4]
    jd = JulianDate(utc=(2014, 1, days))
    pos = earth.at(jd).position.au
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
and other scientific Python tools::

    x, y, z = pos    # four values each
    plot(x, y)

If you instead slice along the second axis,
then you can retrieve an individual position for a particular date —
and the first position is exactly what was returned above
when we computed the January 1st position by itself:

.. testcode::

    print(pos[:,0])

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

Finally, converting an array Julian date back into a calendar tuple
results in the year, month, and all of the other values
being as deep as the array itself:

.. testcode::

    from pprint import pprint
    pprint(jd.utc)

.. testoutput::

    array([[ 2014.,  2014.,  2014.,  2014.],
           [    1.,     1.,     1.,     1.],
           [    1.,     2.,     3.,     4.],
           [    0.,     0.,     0.,     0.],
           [    0.,     0.,     0.,     0.],
           [    0.,     0.,     0.,     0.]])

Again, simply slice across the second dimension of the array
to pull a particular calendar tuple out of the larger result:

.. testcode::

    print(jd.utc[:,2])

.. testoutput::

    [ 2014.     1.     3.     0.     0.     0.]

.. _tai-tt-tdb:

Uniform time scales: TAI, TT, and TDB
=====================================

Date arithmetic becomes very simple
as we leave UTC behind and consider completely uniform time scales.
Days are always twelve hours, hours always 60 minutes,
and minutes always 60 seconds without any variation or dissent.
Such time scales are not appropriate for your morning alarm clock
because they will never be delayed or adjusted
to stay in sync with the slowing rotation of the earth.
But that is what makes them useful for astronomical calculation —
because physics keeps up its dance,
and the stars and planets move in their courses,
whether humanity pauses to observe a UTC leap second or not.

Because they make every day the same length,
uniform time scales can express dates
as a simple floating-point count of days elapsed.
To make all historical dates come out as positive numbers,
astronomers traditionally use Julian days
that start counting at B.C. 4713 January 1 in the old Julian calendar —
the same date as B.C. 4714 November 24 in our Gregorian calendar.
The count starts at noon
following a tradition going back to the Greeks and Ptolemy,
since the sun’s transit is an observable event
but the moment of midnight is not.

So twelve noon was the moment of Julian date zero:

.. testcode::

    bc_4714 = -4713
    print(JulianDate(tt=(bc_4714, 11, 24, 12)).tt)

.. testoutput::

    0.0

Did you notice how negative years work?
People still counted by starting at one, not zero,
when the scholar Dionysius Exiguus created the eras BC and AD
in around the year AD 500.
So his scheme has 1 BC followed immediately by 1 AD without a break.
To avoid an off-by-one error,
astronomers usually ignore BC and count backwards through a year zero
and on into negative years.
So negative year *-n* is what might otherwise be called
either “*n+1* BC” or perhaps “*n+1* BCE” in a history textbook.

More than two million days have passed since 4714 BC
so modern dates tend to be rather large numbers:

.. testcode::

    jd = JulianDate(utc=(2014, 1, 1))
    print('TAI = %r' % jd.tai)
    print('TT  = %r' % jd.tt)
    print('TDB = %r' % jd.tdb)

.. testoutput::

    TAI = 2456658.5004050927
    TT  = 2456658.5007775929
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
that we think an atomic clock would run at the Solar System barycenter,
where it would be unaffected by the Earth’s motion.
The acceleration that Earth experiences in its orbit —
sometimes speeding up, sometimes slowing down —
varies the rate at which our atomic clocks
seem to run to an outside observer,
as predicted by Einstein’s theory of General Relativity.
So physical simulations of the Solar System tend to use TDB,
which is continuous with the T\ :sub:`eph` time scale
traditionally used for Solar System and spacecraft simulations
at the Jet Propulsion Laboratory.

UT1 and ΔT
==========

Finally, UT1 is the least uniform time scale of all
because its clock cannot be housed in a laboratory,
nor is its rate established by any human convention.
It is, rather, the clock
whose “hand” is the rotation of the Earth itself!

The UT1 time is derived from the direction
that the Earth happens to be pointing at any given moment.
And the Earth is a young world
with a still-molten iron core, and viscous mantle,
and continents that rise and fall
as each passing ice age weighs down upon them and then melts away.
We think that we can predict, with high accuracy,
where the planets will be in their orbits
thousands of years from now.
But to predict the fluid dynamics of an elastic rotating ellipsoid
is, at the moment, beyond us.
We cannot run a simulation or formula
to declare leap seconds and keep UTC close to UT1.
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
is provided as ΔT,
the evolving difference between TT and UT1
that dropped below zero in 1871 but then rose past it in 1902
and now stands at more than +67.2 seconds.

The task of governing leap seconds can be stated, then,
as the task of keeping the difference between TT and UTC
close to the wild natural value ΔT.
The standards bodies promise, in fact,
that the difference between these two artificial time scales
will always be within 0.9 seconds of the observed ΔT value.

In calculations that do not involve Earth’s rotation,
ΔT never arises.
The positions of planets,
the distance to the Moon,
and the movement of a comet or asteroid
all ignore ΔT completely.
When, then, does ΔT come into play?

* ΔT is used when you specify your geographic location
  as a :class:`Topos` and Skyfield needs to determine
  where in space that location is facing at a given date and time.

* ΔT is needed to determine directions
  like “up,” “north,” and “east” when you want Skyfield
  to compute the altitude and azimuth of an object
  in your local sky.

* ΔT determines the Earth orientation for Skyfield
  when an Earth satellite position generated from TLE elements
  gets translated into a full Solar System position.

Soon, Skyfield will include a full table of historical ΔT values
along with guidelines about using them in calculations.
For the moment,
if you need the above calculations to run at very high precision,
you can look up ΔT in a table
and provide it to your Julian date manually:

.. testcode::

    JulianDate(utc=(2014, 1, 1), delta_t=67.2810)

.. _date-cache:

The Julian date object as cache
===============================

When you create a :class:`JulianDate`
it goes ahead and computes its ``tt`` Terrestrial Time attribute
starting from whatever time argument you provide.
If you provide the ``utc`` parameter, for example,
then the date first computes and sets ``tai`` followed by ``tt``.
Each of the other time attributes only gets computed once,
the first time you access it.

The general rule is that attributes are only computed once,
and can be accessed again and again for free,
while methods never cache their results —
think of the ``()`` parentheses after a method name
as your reminder that “this will do a fresh computation every time.”

In addition to time scales,
there are several more functions of time
that live on Julian date objects
since they are often needed repeatedly during a calculation.

``gmst``
    Greenwich Mean Sidereal Time.

``gast``
    Greenwich Apparent Sidereal Time.

``P``
    The precession matrix **P** for rotating an *x,y,z* vector
    to the true equator and equinox — the “epoch” — of this Julian date.

``N``
    The even more expensive nutation matrix **N**
    for rotating an *x,y,z* vector to this epoch of this Julian date.

``M``
    The product **NPB** that performs the complete rotation
    between a vector in the ICRS
    and a vector in the dynamical reference system of this Julian date,
    where **B** is the frame tie between the two systems.

``MT``, ``NT``, ``PT``
    The three matrices **M**\ :sup:`T`, **N**\ :sup:`T`,
    and **P**\ :sup:`T`
    that are the transposes of the three previous matrices,
    and that rotate back the other direction
    from the dynamical reference system back to the ICRS frame.

You will typically never have to access these matrices yourself,
as they are used automatically by the :meth:`Position.radec()`
method when you use its  ``epoch=`` parameter
to ask for a right ascension and declination
in the dynamical reference system.

.. _matplotlib: http://matplotlib.org/
.. _pytz: http://pytz.sourceforge.net/
.. _requirements.txt: https://pip.pypa.io/en/latest/user_guide.html#requirements-files
