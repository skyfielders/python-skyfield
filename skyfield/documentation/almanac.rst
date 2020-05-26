
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
    eph = api.load('de421.bsp')

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
    t, y = almanac.find_discrete(t0, t1, almanac.seasons(eph))

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
    t, y = almanac.find_discrete(t0, t1, almanac.moon_phases(eph))

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

.. _lunar-nodes:

Lunar Nodes
===========

The Moon’s ascending node and descending node are the moments each lunar
month when the Moon crosses the plane of Earth’s orbit and eclipses are
possible.

.. testcode::

    t0 = ts.utc(2020, 4, 22)
    t1 = ts.utc(2020, 5, 22)
    t, y = almanac.find_discrete(t0, t1, almanac.moon_nodes(eph))

    print(t.utc_iso())
    print(y)
    print([almanac.MOON_NODES[yi] for yi in y])

.. testoutput::

    ['2020-04-27T17:54:17Z', '2020-05-10T09:01:42Z']
    [ True False]
    ['ascending', 'descending']

.. _oppositions-conjunctions:

Opposition and Conjunction
==========================

The moment at which a planet is in opposition with the Sun or in
conjunction with the Sun is when their ecliptic longitudes are at 0° or
180° difference.

.. testcode::

    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2021, 1, 1)
    f = almanac.oppositions_conjunctions(eph, eph['mars'])
    t, y = almanac.find_discrete(t0, t1, f)

    print(t.utc_iso())
    print(y)

.. testoutput::

    ['2019-09-02T10:42:14Z', '2020-10-13T23:25:47Z']
    [0 1]

The result ``t`` will be an array of times, and ``y`` will be an array
of integers indicating which half of the sky the body has just entered:
0 means the half of the sky west of the Sun along the ecliptic, and 1
means the half of the sky east of the Sun.  This means different things
for different bodies:

* For the outer planets Mars, Jupiter, Saturn, Uranus, and all other
  bodies out beyond our orbit, 0 means the moment of conjunction with
  the Sun and 1 means the moment of opposition.

* Because the Moon moves eastward across our sky relative to the Sun,
  not westward, the output is reversed compared to the outer planets: 0
  means the moment of opposition or Full Moon, while 1 means the moment
  of conjunction or New Moon.

* The inner planets Mercury and Venus only ever experience conjunctions
  with the Sun from our point of view, never oppositions, with 0
  indicating an inferior conjunction and 1 a superior conjunction.

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
    t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(eph, bluffton))

    print(t.utc_iso())
    print(y)

.. testoutput::

    ['2018-09-12T11:13:13Z', '2018-09-12T23:49:38Z']
    [ True False]

The result ``t`` will be an array of times, and ``y`` will be ``True``
if the sun rises at the corresponding time and ``False`` if it sets.

If you need to provide your own custom value for refraction, adjust the
estimate of the Sun’s radius, or account for a vantage point above the
Earth’s surface, see :ref:`risings-and-settings` to learn about the more
versatile :func:`~skyfield.almanac.risings_and_settings()` routine.

Note that a location near one of the poles during polar summer or polar
winter will not experience sunrise and sunset.  To learn whether the sun
is up or down, call the sunrise-sunset function at the time that
interests you, and the return value will indicate whether the sun is up.

.. testcode::

    far_north = api.Topos('89 N', '80 W')
    f = almanac.sunrise_sunset(eph, far_north)
    t, y = almanac.find_discrete(t0, t1, f)

    print(t.utc_iso())  # Empty list: no sunrise or sunset
    print(f(t0))        # But we can ask if the sun is up

    print('polar day' if f(t0) else 'polar night')

.. testoutput::

    []
    True
    polar day

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
    f = almanac.dark_twilight_day(eph, bluffton)
    t, y = almanac.find_discrete(t0, t1, f)

    for ti, yi in zip(t, y):
        print(yi, ti.utc_iso(), ' Start of', almanac.TWILIGHTS[yi])

.. testoutput::

    1 2019-11-08T10:40:19Z  Start of Astronomical twilight
    2 2019-11-08T11:12:31Z  Start of Nautical twilight
    3 2019-11-08T11:45:18Z  Start of Civil twilight
    4 2019-11-08T12:14:15Z  Start of Day
    3 2019-11-08T22:23:52Z  Start of Civil twilight
    2 2019-11-08T22:52:48Z  Start of Nautical twilight
    1 2019-11-08T23:25:34Z  Start of Astronomical twilight
    0 2019-11-08T23:57:44Z  Start of Night

.. _risings-and-settings:

Risings and Settings
====================

Skyfield can compute when a given body rises and sets.
The routine is designed for bodies at the Moon’s distance or farther,
that tend to rise and set about once a day.
But it might be caught off guard
if you pass it an Earth satellite
that rises several times a day;
for that case, see :ref:`satellite-rising-and-setting`.

Rising and setting predictions can be generated
using the :func:`~skyfield.almanac.risings_and_settings()` routine:

.. testcode::

    t0 = ts.utc(2020, 2, 1)
    t1 = ts.utc(2020, 2, 2)
    f = almanac.risings_and_settings(eph, eph['Mars'], bluffton)
    t, y = almanac.find_discrete(t0, t1, f)

    for ti, yi in zip(t, y):
        print(ti.utc_iso(), 'Rise' if yi else 'Set')

.. testoutput::

    2020-02-01T09:29:16Z Rise
    2020-02-01T18:42:57Z Set

As with sunrise and sunset above, ``True`` means the moment of rising
and ``False`` means the moment of setting.

The routine also offers some optional parameters,
whose several uses are covered in the following sections.

Computing your own refraction angle
-----------------------------------

Instead of accepting the standard estimate of 34 arcminutes
for the angle by which refraction will raise the image
of a body at the horizon,
you can compute atmospheric refraction yourself
and supply the resulting angle to ``horizon_degrees``.
Note that the value passed should be a small negative angle.
In this example it makes a 3 second difference
in both the rising and setting time:

.. testcode::

    from skyfield.earthlib import refraction

    r = refraction(0.0, temperature_C=15.0, pressure_mbar=1030.0)
    print('Arcminutes refraction for body seen at horizon: %.2f\n' % (r * 60.0))

    f = almanac.risings_and_settings(eph, eph['Mars'], bluffton, horizon_degrees=-r)
    t, y = almanac.find_discrete(t0, t1, f)

    for ti, yi in zip(t, y):
        print(ti.utc_iso(), 'Rise' if yi else 'Set')

.. testoutput::

    Arcminutes refraction for body seen at horizon: 34.53

    2020-02-01T09:29:13Z Rise
    2020-02-01T18:43:00Z Set

Adjusting for apparent radius
-----------------------------

Planets and especially the Sun and Moon have an appreciable radius,
and we usually consider the moment of sunrise
to be the moment when its bright limb crests the horizon —
not the later moment when its center finally rises into view.
Set the parameter ``radius_degrees`` to the body’s apparent radius
to generate an earlier rising and later setting;
the value ``0.25``, for example,
would be a rough estimate for the Sun or Moon.

The difference in rising time can be a minute or more:

.. testcode::

    f = almanac.risings_and_settings(eph, eph['Sun'], bluffton, radius_degrees=0.25)
    t, y = almanac.find_discrete(t0, t1, f)
    print(t[0].utc_iso(' '), 'Limb of the Sun crests the horizon')

    f = almanac.risings_and_settings(eph, eph['Sun'], bluffton)
    t, y = almanac.find_discrete(t0, t1, f)
    print(t[0].utc_iso(' '), 'Center of the Sun reaches the horizon')

.. testoutput::

    2020-02-01 12:46:28Z Limb of the Sun crests the horizon
    2020-02-01 12:47:53Z Center of the Sun reaches the horizon

Elevated vantage points
-----------------------

Rising and setting predictions usually assume a flat local horizon
that does not vary with elevation.
Yes, Denver is the Mile High City,
but it sees the sun rise against a local horizon that’s also a mile high.
Since the city’s high altitude
is matched by the high altitude of the terrain around it,
the horizon winds up in the same place it would be for a city at sea level.

But sometimes you need to account not only for local elevation,
but for *altitude* above the surrounding terrain.
Some observatories, for example, are located on mountaintops
that are much higher than the elevation of the terrain
that forms their horizon.
And Earth satellites can be hundreds of kilometers
above the surface of the Earth that produces their sunrises and sunsets.

You can account for high altitude above the horizon’s terrain
by setting an artificially negative value for ``horizon_degrees``.
If we consider the Earth to be approximately a sphere,
then we can use a bit of trigonometry
to estimate the position of the horizon for an observer at altitude:

.. testcode::

    from numpy import arccos
    from skyfield.units import Angle

    # When does the Sun rise in the ionosphere’s F-layer, 300km up?
    altitude_m = 300e3

    earth_radius_m = 6378136.6
    side_over_hypotenuse = earth_radius_m / (earth_radius_m + altitude_m)
    h = Angle(radians = -arccos(side_over_hypotenuse))
    print('The horizon from 300km up is at %.2f degrees' % h.degrees)

    f = almanac.risings_and_settings(
        eph, eph['Sun'], bluffton, horizon_degrees=h.degrees,
        radius_degrees=0.25,
    )
    t, y = almanac.find_discrete(t0, t1, f)
    print(t[0].utc_iso(' '), 'Limb of the Sun crests the horizon')

.. testoutput::

    The horizon from 300km up is at -17.24 degrees
    2020-02-01 00:22:42Z Limb of the Sun crests the horizon

When writing code for this situation,
we need to be very careful to keep straight
the two different meanings of *altitude*.

1. The *altitude above sea level* is a linear distance measured in meters
   between the ground and the location at which
   we want to compute rises and settings.

2. The *altitude of the horizon* names a quite different measure.
   It’s an angle measured in degrees
   that is one of the two angles
   of the altitude-azimuth (“altazimuth”) system
   oriented around an observer on a planet’s surface.
   While azimuth measures horizontally around the horizon
   from north through east, south, and west,
   the altitude angle measures up towards the zenith (positive)
   and down towards the nadir (negative).
   The altitude is zero all along the great circle between zenith and nadir.

The problem of an elevated observer
unfortunately involves both kinds of altitude at the same time:
for each extra meter of “altitude” above the ground,
there is a slight additional depression in the angular “altitude”
of the horizon on the altazimuth globe.

When a right ascension and declination rises and sets
-----------------------------------------------------

If you are interested in finding the times
when a fixed point in the sky rises and sets,
simply create a star object with the coordinates
of the position you are interested in
(see :doc:`stars`).
Here, for example, are rising and setting times for the Galactic Center:

.. testcode::

    galactic_center = api.Star(ra_hours=(17, 45, 40.04),
                               dec_degrees=(-29, 0, 28.1))

    f = almanac.risings_and_settings(eph, galactic_center, bluffton)
    t, y = almanac.find_discrete(t0, t1, f)

    for ti, yi in zip(t, y):
        verb = 'rises above' if yi else 'sets below'
        print(ti.utc_iso(' '), '- Galactic Center', verb, 'the horizon')

.. testoutput::

    2020-02-01 10:29:00Z - Galactic Center rises above the horizon
    2020-02-01 18:45:46Z - Galactic Center sets below the horizon

Solar terms
===========

The solar terms are widely used in East Asian calendars.

.. testcode::

    from skyfield import almanac_east_asia as almanac_ea

    t0 = ts.utc(2019, 12, 1)
    t1 = ts.utc(2019, 12, 31)
    t, tm = almanac.find_discrete(t0, t1, almanac_ea.solar_terms(eph))

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
