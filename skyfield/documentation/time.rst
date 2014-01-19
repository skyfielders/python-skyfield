
.. currentmodule:: skyfield.api

================
 Dates and Time
================

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

* ``jd.utc_…`` — Coordinated Universal Time (“Greenwich Time”)
* ``jd.tai`` — International Atomic Time
* ``jd.tt`` — Terrestrial Time
* ``jd.tdb`` — Barycentric Dynamical Time (the JPL’s “T\ :sub:`eph`”)
* ``jd.ut1`` — Universal Time

To build a Julian date object,
simply use one of these time scales
as a keyword argument to specify the moment you want to represent:

.. testcode::

    from skyfield.api import JulianDate

    jd = JulianDate(utc=(2014, 1, 18))

For all of the techniques that you can use to build a Julian date,
read the next two sections on :ref:`utc-and-timezone`
and on :ref:`building-dates`.

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

    from skyfield.api import earth

    jd = JulianDate(utc=(2014, 1, 1))
    print earth(jd).position.AU

    print earth(utc=(2014, 1, 1)).position.AU

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

If you are comfortable dealing with UTC,
you can simply set and retrieve it manually.
A Python tuple is the most convenient way
to represent the year, month, and day of a calendar date,
and you have the option of providing hours, minutes, and seconds
as well:

.. testcode::

    from skyfield.api import utc

    # Four ways to specify 2014 January 18 01:35:37.5

    jd  = JulianDate(utc=(2014, 1, 18.06640625))
    jd2 = JulianDate(utc=(2014, 1, 18, 1.59375))
    jd3 = JulianDate(utc=(2014, 1, 18, 1, 35.625))
    jd4 = JulianDate(utc=(2014, 1, 18, 1, 35, 37.5))

    assert jd == jd2 == jd3 == jd4

    # Retrieving UTC as either a tuple or string

    print jd.utc
    print jd.utc_iso()
    print jd.utc_jpl()
    print jd.utc_strftime('Date %Y-%m-%d and time %H:%M:%S')

.. testoutput::

    (2014, 1, 18, 1.0, 35.0, 37.5)
    2014-01-18T01:35:38Z
    A.D. 2014-Feb-18 01:35:37.5000 UT
    Date 2014-01-18 and time 01:35:38

And by scraping together the minimal support for UTC
that exists in the Python Standard Library,
Skyfield is able to return the current time as a Julian date object
(assuming that your operating system clock is correct
and is configured with the correct time zone):

.. testsetup::

    from skyfield import api
    def now():
        """Return a constant "now"."""
        return api.JulianDate(utc=(2014, 1, 18, 23, 10, 9))
    api.now = now

.. testcode::

    from skyfield.api import now

    jd = now()
    print jd.utc_jpl()

.. testoutput::

    A.D. 2014-Feb-18 23:10:09.0000 UT

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
    from pytz import timezone, utc

    eastern = timezone('US/Eastern')

    # Converting US Eastern Time to a Julian date.

    d = datetime(2014, 1, 16, 1, 32, 9)
    e = eastern.localize(d)
    jd = JulianDate(utc=e)

    # Building a datetime in UTC

    dt, leap_second = jd.utc_datetime()
    print 'UTC:', dt

    # Converting back to an Eastern Time datetime.

    dt, leap_second = jd.astimezone(eastern)
    print 'EST:', dt

.. testoutput::

    UTC: 2014-01-16 06:32:09+00:00
    EST: 2014-01-16 01:32:09-05:00

As we would expect,
1:32 AM in the Eastern time zone in January
is 6:32 AM local time in Greenwich, England,
five hours to the east across the Atlantic.
(The extra return value ``leap_second``
is explained in the next section, :ref:`leap-seconds`.)

Note that the :meth:`~JulianDate.astimezone()` method
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
a day of 24 hours, and an hour of 60 minutes,
the only way to stay within the rules
while keeping UTC synchronized with the Earth
is to occasionally add an extra leap second
to one of the year’s minutes.

The `International Earth Rotation Service <http://hpiers.obspm.fr/>`_
currently restricts itself to adding a leap second
into the last minute of June or the last minute of December.
When a leap second is inserted,
its minute counts 61 seconds numbered 00–60
instead of staying within the usual range 00–59.
The most recent leap second was in June 2012:

.. testcode::

    five_seconds = range(58, 58 + 5)
    tup = (2012, 6, 30, 23, 59, five_seconds)
    jd = JulianDate(utc=tup)

    for string in jd.utc_jpl():
        print string

.. testoutput::

    A.D. 2012-Jul-30 23:59:58.0000 UT
    A.D. 2012-Jul-30 23:59:59.0000 UT
    A.D. 2012-Jul-30 23:59:60.0000 UT
    A.D. 2012-Aug-01 00:00:00.0000 UT
    A.D. 2012-Aug-01 00:00:01.0000 UT

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
from use as a general way to represent time,
because it refuses to accept seconds greater than 59:

.. testcode::

    datetime(2012, 6, 30, 19, 59, 60)

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: second must be in 0..59

That is why operations that return a ``datetime``
always include a second return value, ``leap_second``,
which is usually ``0`` but jumps to ``1``
when Skyfield is forced to represent a leap second
as a duplicate 23:59:59, as is the case here:

.. testcode::

    dt_list, leap_seconds = jd.astimezone(eastern)

    for dt, leap_second in zip(dt_list, leap_seconds):
        print str(dt) + (' +1s' if leap_second else '')

.. testoutput::

    2012-06-30 19:59:58-04:00
    2012-06-30 19:59:59-04:00
    2012-06-30 19:59:59-04:00 +1s
    2012-06-30 20:00:00-04:00
    2012-06-30 20:00:01-04:00

Using tuples to represent UTC times is more elegant,
because leap seconds can be represented accurately.
If your application cannot avoid using ``datetime`` objects,
then you will have to decide
whether to simply ignore the ``leap_second`` value
or to somehow output the leap second information.

Note that it is always preferred to use the ``leap_second`` name
in your code even when you do nothing with the value,
to document that you have made a conscious choice
to ignore this particular complication:

.. testcode::

    # Bad: your code fails to document what it is ignoring.

    dt_list = jd.astimezone(eastern)[0]

    # Good: the meaning of the second value is explicit,
    # even if you make no use of it.

    dt_list, leap_seconds = jd.astimezone(eastern)

.. _date-arrays:

Date arrays
===========

To compute a position or coordinate
over a whole range of dates,
take advantage of Skyfield’s support for NumPy arrays.
Building separate date objects is a slow process:

.. testcode::

    # Building separate dates is slow and awkward

    for day in 1, 2, 3, 4:
        jd = JulianDate(utc=(2014, 1, day))
        print earth(jd).position.AU

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]
    [-0.19179872  0.88265548  0.38254134]
    [-0.20891924  0.87936337  0.38111391]
    [-0.22597338  0.87579547  0.37956709]

If you instead construct a single :class:`JulianDate`
that holds a whole array of dates,
then explicit loops disappear from your code
and Skyfield can perform efficient vector computations
over the whole range of dates at once:

.. testcode::

    # Date arrays are concise and efficient

    days = [1, 2, 3, 4]
    jd = JulianDate(utc=(2014, 1, days))
    pos = earth(jd).position.AU
    print pos

.. testoutput::

    [[-0.17461758 -0.19179872 -0.20891924 -0.22597338]
     [ 0.88567056  0.88265548  0.87936337  0.87579547]
     [ 0.38384886  0.38254134  0.38111391  0.37956709]]

Note the shape of the resulting NumPy array.
If you unpack this array using three names,
then you get three four-element arrays
whose four values correspond
to the four dates in the :class:`JulianDate`.
These arrays are ready to be submitted to `matplotlib`_
and other scientific Python tools::

    x, y, z = pos
    plot(x, y)

If you instead slice along the second axis,
then you can retrieve an individual position for a particular date —
for example, the position corresponding to the first date
lives at index zero:

.. testcode::

    print pos[:,0]

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

If you want to check,
you can scroll up and verify that this is the same coordinate
that we generated for January 1st at the top of this document
when we were using a single date instead of an array.

.. _utc-and-timezone:

Building UTC arrays
===================



Other timescales: TAI, TT, and TDB
==================================

.. testcode::

    jd = JulianDate(utc=(2014, 1, 1))
    print 'TAI = %r' % jd.tai
    print 'TT  = %r' % jd.tt
    print 'TDB = %r' % jd.tdb
    print 'UT1 = %r' % jd.ut1

.. testoutput::

    TAI = 2456658.5004050927
    TT  = 2456658.5007775929
    TDB = 2456658.500777592
    UT1 = 2456658.5007775929

UT1 and ΔT
==========

.. _date-cache:

The Julian date object as cache
===============================


.. _matplotlib: http://matplotlib.org/
.. _pytz: http://pytz.sourceforge.net/
.. _requirements.txt: http://www.pip-installer.org/en/latest/cookbook.html
