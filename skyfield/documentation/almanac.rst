
=====================
 Almanac Computation
=====================

The highest-level routines in Skyfield let you search back and forward
through time for the exact moments when the Earth, Sun, and Moon are in
special configurations.

They all require you to start by loading up a timescale object and also
an ephemeris file that provides positions from the planets:

.. testcode::

    from skyfield import api

    ts = api.load.timescale()
    e = api.load('de421.bsp')

Then, load the “almanac” module.

.. testcode::

    from skyfield import almanac

Note that almanac computation can be slow and expensive.  To determine
the moment of sunrise, for example, Skyfield has to search back and
forth through time asking for the altitude of the Sun over and over
until it finally works out the moment at which it crests the horizon.

Rounding time to the nearest minute
===================================

If you compare almanac results to official sources like the `United
States Naval Observatory <http://aa.usno.navy.mil/data/index.php>`_, the
printed time will often differ because the Naval Observatory results are
rounded to the nearest minute — any time with ``:30`` or more seconds at
the end gets named as the next minute.

If you try to display a date that needs to be rounded to the nearest
minute by simply stopping at ``%M`` and leaving off the ``%S`` seconds,
the output will be one minute too early.  For example, the Naval
Observatory would round ``14:59`` up to ``:15`` in the following date.

.. testcode::

    t = ts.utc(2018, 9, 10, 5, 14, 59)
    dt = t.utc_datetime()
    print(dt.strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    2018-09-10 05:14

To do the same rounding yourself, simply add 30 seconds to the time
before truncating the seconds.

.. testcode::

    from datetime import timedelta

    def nearest_minute(dt):
        return (dt + timedelta(seconds=30)).replace(second=0, microsecond=0)

    dt = nearest_minute(t.utc_datetime())
    print(dt.strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    2018-09-10 05:15

The results should then agree with the tables produced by the USNO.

The Seasons
===========

Create a start time and an end time to ask for all of the equinoxes and
solstices that fall in between.

.. testcode::

    t0 = ts.utc(2018, 1, 1)
    t1 = ts.utc(2018, 12, 31)
    t, y = almanac.find_discrete(t0, t1, almanac.seasons(e))

    for yi, ti in zip(y, t):
        print(yi, almanac.SEASON_EVENTS[yi], ti.utc_iso(' '))

.. testoutput::

    0 Vernal Equinox 2018-03-20 16:15:27Z
    1 Summer Solstice 2018-06-21 10:07:18Z
    2 Autumnal Equinox 2018-09-23 01:54:06Z
    3 Winter Solstice 2018-12-21 22:22:44Z

The result ``t`` will be an array of times, and ``y`` will be ``0``
through ``3`` for the Vernal Equinox through the Winter Solstice.


Sunrise and Sunset
==================

Because sunrise and sunset differ depending on your location on the
Earth’s surface, you first need to create a Topos object describing your
geographic location.

.. testcode::

    bluffton = api.Topos('40.8939 N', '83.8917 W')

Then you can create a start time and an end time and ask for all of the
sunrises and sunsets in between.

.. testcode::

    t0 = ts.utc(2018, 9, 12, 4)
    t1 = ts.utc(2018, 9, 13, 4)
    t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(e, bluffton))

    print(t.utc_iso())
    print(y)

.. testoutput::

    ['2018-09-12T11:13:13Z', '2018-09-12T23:49:38Z']
    [ True False]

The result ``t`` will be an array of times, and ``y`` will be ``True``
if the sun rises at the corresponding time and ``False`` if it sets.

Phases of the Moon
==================

The phases of the Moon are the same for everyone on Earth, so no Topos
is necessary but only an ephemeris object.

.. testcode::

    t0 = ts.utc(2018, 9, 1)
    t1 = ts.utc(2018, 9, 10)
    t, y = almanac.find_discrete(t0, t1, almanac.moon_phases(e))

    print(t.utc_iso())
    print(y)
    print([almanac.MOON_PHASES[yi] for yi in y])

.. testoutput::

    ['2018-09-03T02:37:24Z', '2018-09-09T18:01:28Z']
    [3 0]
    ['Last Quarter', 'New Moon']

The result ``t`` will be an array of times, and ``y`` will be a
corresponding array of Moon phases with 0 for New Moon and 3 for Last
Quarter.  You can use the array ``MOON_PHASES`` to retrieve names for
each phase.
