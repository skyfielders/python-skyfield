
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

* ``jd.utc`` — Coordinated Universal Time (“Greenwich Time”)
* ``jd.tai`` — International Atomic Time
* ``jd.tt`` — Terrestrial Time
* ``jd.tdb`` — Barycentric Dynamical Time (the JPL’s “T\ :sub:`eph`”)
* ``jd.ut1`` — Universal Time

To build a Julian date object,
simply use one of these time scales
as a keyword argument to the constructor
to specify the moment you want to represent:

.. testcode::

    from skyfield.api import JulianDate

    jd = JulianDate(utc=(2014, 1, 1))

For all of the techniques that you can use to build a Julian date,
read the next two sections on :ref:`your-timezone`
and on :ref:`building-dates`.

There are two ways to provide a date argument
to a Skyfield routine that needs one:
you can supply a :class:`JulianDate`
that you have gone ahead and built yourself;
or you can simply call the routine with the same keyword argument
that you would have used to build the Julian date,
and Skyfield will construct the date internally.
So these two calculations come out exactly the same:

.. testcode::

    from skyfield.api import earth

    print earth(jd).position.AU
    print earth(utc=(2014, 1, 1)).position.AU

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]
    [-0.17461758  0.88567056  0.38384886]

Be sure to build the Julian date yourself
if you have the opportunity to use it over and over again
in several different calculations.
Not only will you avoid recomputing the time scales
again for each calculation,
but there are some even more expensive functions of time
that Julian dates take responsibility for caching
that you can read about below in the section on :ref:`date-cache`.

.. _building-dates:

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

.. _your-timezone:

Your timezone
=============

The :class:`JulianDate` object in this document’s first example
represented midnight on 2014 January 1
in Coordinated Universal Time (UTC)
which is the world clock — once known as “Greenwich Mean Time” —
from which all other time zones are mere offsets.

The Python Standard Library provides a
`datetime.utcnow() <http://docs.python.org/2/library/datetime.html#datetime.datetime.utcnow>`_
constructor that returns the current UTC date and time,
if your operating system time zone is correctly configured.
Skyfield’s :class:`JulianDate` knows how to parse a ``datetime``
if you supply it as the UTC argument:

.. testcode::

    from datetime import datetime
    jd = JulianDate(utc=datetime.utcnow())

But to perform actual conversion to and from your own timezone,
you will need to install the `pytz <http://pytz.sourceforge.net/>`_
package, either by listing it in the dependencies of your package,
adding it to your project’s `requirements.txt`_ file,
or simply installing it manually::

    pip install pytz

Building dates is simple:
instantiate a normal Python ``datetime`` object,
and then pass it to the ``localize()`` method
of your time zone.
You can then pass the date safely to Skyfield
and the result will come out correctly as a UTC time:

.. testcode::

    from pytz import timezone, utc
    eastern = timezone('US/Eastern')

    d = datetime(2014, 1, 16, 1, 32, 9)
    e = eastern.localize(d)
    jd = JulianDate(utc=e)
    print jd.utc_datetime()

.. testoutput::

    2014-01-16 06:32:09+00:00

As we would expect,
1:32 AM in the Eastern time zone in January
is 6:32 AM — just an hour before dawn — in Greenwich, England,
five hours to the east across the Atlantic.

To convert back to your own timezone,
use the :meth:`JulianDate.astimezone()` method.
When Skyfield sees that your timezone is a ``pytz`` timezone,
it will even call its ``normalize()`` method for you
to make sure that daylight savings time is handled correctly!

.. testcode::

    print jd.astimezone(eastern)

.. testoutput::

    2014-01-16 01:32:09-05:00

And because :class:`JulianDate` fully supports arrays,
all of these techniques work just fine
with whole lists of ``datetime`` objects as inputs and outputs:

.. testcode::

    hours = range(4)
    dates = [eastern.localize(datetime(2014, 1, 16, h))
             for h in hours]

    jd = JulianDate(utc=dates)

    from pprint import pprint
    pprint(jd.utc_datetime())

.. testoutput::

    [datetime.datetime(2014, 1, 16, 5, 0, tzinfo=<UTC>),
     datetime.datetime(2014, 1, 16, 6, 0, tzinfo=<UTC>),
     datetime.datetime(2014, 1, 16, 7, 0, tzinfo=<UTC>),
     datetime.datetime(2014, 1, 16, 8, 0, tzinfo=<UTC>)]

.. testcode::

    pprint(jd.astimezone(eastern))

.. testoutput::

    [datetime.datetime(2014, 1, 16, 0, 0, tzinfo=...EST...),
     datetime.datetime(2014, 1, 16, 1, 0, tzinfo=...EST...),
     datetime.datetime(2014, 1, 16, 2, 0, tzinfo=...EST...),
     datetime.datetime(2014, 1, 16, 3, 0, tzinfo=...EST...)]

Building Python lists of ``datetime`` objects
is neither as elegant nor as efficient
as working directly in UTC with NumPy arrays,
which we will explore in more detail the next section.
But if you need to operate primarily in your own timezone,
then it is the only option.

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
.. _requirements.txt: http://www.pip-installer.org/en/latest/cookbook.html
