
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

If you or some of your users live in the Southern Hemisphere,
you can use the ``SEASON_EVENTS_NEUTRAL`` array.
Instead of naming specific seasons,
it names the equinoxes and solstices by the month in which they occur —
so the ``March Equinox``, for example, is followed by the ``June Solstice``.

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

Opposition and Conjunction
==========================

The moment at which a planet is in opposition with the Sun or in
conjunction with the Sun is when their ecliptic longitudes are at 0° or
180° difference.

.. testcode::

    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2021, 1, 1)
    f = almanac.oppositions_conjunctions(e, e['mars'])
    t, y = almanac.find_discrete(t0, t1, f)

    print(t.utc_iso())
    print(y)
    print([almanac.CONJUNCTIONS[yi] for yi in y])

.. testoutput::

    ['2019-09-02T10:42:14Z', '2020-10-13T23:25:47Z']
    [0 1]
    ['conjunction', 'opposition']

The result ``t`` will be an array of times, and ``y`` will be an array
of integers where 0 means a conjunction and 1 means an opposition.  You
can use the array ``CONJUNCTIONS`` to retrieve names for each value.

Sunrise and Sunset
==================

Because sunrise and sunset differ depending on your location on the
Earth’s surface, you first need to create a Topos object describing your
geographic location.

.. testcode::

    bluffton = api.Topos('40.8939 N', '83.8917 W')

Then you can create a start time and an end time and ask for all of the
sunrises and sunsets in between.
Skyfield uses the
`official definition of sunrise and sunset
<http://aa.usno.navy.mil/faq/docs/RST_defs.php>`_
from the United States Naval Observatory,
which defines them as the moment when the center — not the limb —
of the sun is 0.8333 degrees below the horizon,
to account for both the average radius of the Sun itself
and for the average refraction of the atmosphere at the horizon.

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

Twilight
========

An expanded version of the sunrise-sunset routine separately codes each
of the phases of twilight using integers:

0. Dark of night.
1. Astronomical twilight.
2. Nautical twilight.
3. Civil twilight.
4. Daytime.

.. testcode::

    t0 = ts.utc(2019, 11, 8, 5)
    t1 = ts.utc(2019, 11, 9, 5)
    t, y = almanac.find_discrete(t0, t1, almanac.dark_twilight_day(e, bluffton))
    for ti, yi in zip(t, y):
        print(yi, ti.utc_iso())

.. testoutput::

    1 2019-11-08T10:40:20Z
    2 2019-11-08T11:12:31Z
    3 2019-11-08T11:45:18Z
    4 2019-11-08T12:14:15Z
    3 2019-11-08T22:23:52Z
    2 2019-11-08T22:52:49Z
    1 2019-11-08T23:25:34Z
    0 2019-11-08T23:57:44Z

Solar terms
===========

The solar terms are widely used in East Asian calendars.

.. testcode::

    from skyfield import almanac_east_asia as almanac_ea

    t0 = ts.utc(2019, 12, 1)
    t1 = ts.utc(2019, 12, 31)
    t, tm = almanac.find_discrete(t0, t1, almanac_ea.solar_terms(e))

    for tmi, ti in zip(tm, t):
        print(tmi, almanac_ea.SOLAR_TERMS_ZHS[tmi], ti.utc_iso(' '))

.. testoutput::

    17 大雪 2019-12-07 10:18:28Z
    18 冬至 2019-12-22 04:19:26Z

The result ``t`` will be an array of times, and ``y`` will be integers
in the range 0–23 which are each the index of a solar term.  Localized
names for the solar terms in different East Asia languages are provided
as ``SOLAR_TERMS_JP`` for Japanese, ``SOLAR_TERMS_VN`` for Vietnamese,
``SOLAR_TERMS_ZHT`` for Traditional Chinese, and (as shown above)
``SOLAR_TERMS_ZHS`` for Simplified Chinese.
