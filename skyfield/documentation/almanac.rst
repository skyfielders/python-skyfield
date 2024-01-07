
=====================
 Almanac Computation
=====================

To help you navigate, here are the topics covered on this page:

.. contents::
   :local:

Here are the imports and objects that will drive the examples below:

.. testcode::

    from skyfield import almanac
    from skyfield.api import N, W, load, wgs84

    ts = load.timescale()
    eph = load('de421.bsp')
    bluffton = wgs84.latlon(40.8939 * N, 83.8917 * W)
    observer = eph['Earth'] + bluffton

Note a distinction:

* ``bluffton`` — this :class:`~skyfield.toposlib.GeographicPosition`
  computes the vector from the Earth’s center
  to a specific geographic location.
  Unless an elevation is specified,
  the vector’s length will be the Earth’s radius,
  which varies from 6,357 km at the pole to 6,378 km at the equator.

* ``observer`` — this computes the much larger vector
  from the center of the Solar System to the geographic position.

  Almanac routines typically need to be passed this full vector,
  not merely a geographic position.
  Depending on the season (Earth is closest to the Sun in early January),
  it will be somewhere close to 1 au in length.

Beware that almanac computation can be expensive.
As they search backwards and forwards through time
for particular circumstances,
almanac routines have to repeatedly recompute
the positions of the Earth and other bodies.

Rounding time to the nearest minute
===================================

If you compare almanac results to official sources like the `United
States Naval Observatory <https://aa.usno.navy.mil/data/index>`_, the
printed time will often differ because the Naval Observatory results are
rounded to the nearest minute — any time with ``:30`` or more seconds at
the end gets rounded up to the next minute.

The Skyfield method :meth:`~skyfield.timelib.Time.utc_strftime()`
performs this rounding automatically if you don’t ask for the seconds to
be displayed:

.. testcode::

    t = ts.utc(2023, 8, 10, 6, 21, 45.9)

    print('Microseconds:', t.utc_strftime('%Y-%m-%d %H:%M:%S.%f'))
    print('Nearest second:', t.utc_strftime('%Y-%m-%d %H:%M:%S'))
    print('Nearest minute:', t.utc_strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    Microseconds: 2023-08-10 06:21:45.900000
    Nearest second: 2023-08-10 06:21:46
    Nearest minute: 2023-08-10 06:22

But beware that the rounding won’t happen automatically if you are doing
your own formatting using Python’s built-in ``datetime`` objects.  For
example, if you stop with ``%M``, then the seconds are simply ignored,
instead of being used for rounding:

.. testcode::

    dt = t.utc_datetime()
    print(dt.strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    2023-08-10 06:21

To fix the problem and round a Python ``datetime`` to the nearest
minute, try manually adding 30 seconds to the time before displaying it:

.. testcode::

    from datetime import timedelta

    def nearest_minute(dt):
        return (dt + timedelta(seconds=30)).replace(second=0, microsecond=0)

    dt = nearest_minute(t.utc_datetime())
    print(dt.strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    2023-08-10 06:22

The results should then agree with the tables produced by the USNO.

.. _risings-and-settings:

Risings and settings
====================

Skyfield can compute when a given body rises and sets
for an observer at the Earth’s surface.
The routine is designed for bodies at the Moon’s distance or farther,
that rise and set about once a day,
so it will be caught off-guard
if you pass it something fast like an Earth satellite;
for that case, see :ref:`satellite-rising-and-setting`.

Sunrise and Sunset
------------------

Skyfield uses the
`official definition of sunrise and sunset
<https://aa.usno.navy.mil/faq/RST_defs>`_
from the United States Naval Observatory,
which defines them as the moment when the center
of the sun is 0.8333° below the horizon,
to account for both the average radius of the Sun
and for the average refraction of the atmosphere at the horizon.
Here’s how to ask for the sunrises between a given start and end time:

.. testcode::

    observer = eph['Earth'] + bluffton
    sun = eph['Sun']
    t0 = ts.utc(2018, 9, 12, 4)
    t1 = ts.utc(2018, 9, 14, 4)

    t, y = almanac.find_risings(observer, sun, t0, t1)
    print(t.utc_iso(' '))
    print(y)

.. testoutput::

    ['2018-09-12 11:13:12Z', '2018-09-13 11:14:12Z']
    [ True  True]

And here’s how to ask for the sunsets:

.. testcode::

    t, y = almanac.find_settings(observer, sun, t0, t1)
    print(t.utc_iso(' '))
    print(y)

.. testoutput::

    ['2018-09-12 23:49:38Z', '2018-09-13 23:47:56Z']
    [ True  True]

Normally every value in the second array will be ``True``,
indicating that a rising or setting was successfully detected.
But locations north of the Arctic Circle or south of the Antarctic Circle
can experience 24-hour summer days during which the sun never sets,
and suffer winter days during which the sun never rises.
On days when the Sun never reaches the horizon,
the second array will have the value ``False``,
and Skyfield will return the moment of transit instead.
For example:

.. testcode::

    harra_sweden = wgs84.latlon(+67.4066, +20.0997)
    harra_observer = eph['Earth'] + harra_sweden
    sun = eph['Sun']

    t0 = ts.utc(2022, 12, 18)
    t1 = ts.utc(2022, 12, 26)
    t, y = almanac.find_risings(harra_observer, sun, t0, t1)

    alt, az, dist = harra_observer.at(t).observe(sun).apparent().altaz()

    for ti, yi, alti in zip(t.utc_iso(' '), y, alt.degrees):
        print('{} {:5} {:.4f}'.format(ti, str(yi), alti))

.. testoutput::

    2022-12-18 10:22:54Z True  -0.8333
    2022-12-19 10:29:21Z True  -0.8333
    2022-12-20 10:37:06Z False -0.8387
    2022-12-21 10:37:36Z False -0.8464
    2022-12-22 10:38:06Z False -0.8461
    2022-12-23 10:38:36Z False -0.8380
    2022-12-24 10:31:28Z True  -0.8333
    2022-12-25 10:26:08Z True  -0.8333

This output shows that right around the winter solstice,
there are four days on which the Sun never quite reaches the horizon,
but is at least a few fractions of a degree below the altitude of -0.8333°
that would qualify for the USNO definition of sunrise.
So Skyfield instead returns the moment when the Sun is closest to the horizon,
with the accompanying value ``False``.



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

As with sunrise and sunset above,
``1`` means the moment of rising and ``0`` means the moment of setting.

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

    2020-02-01 12:46:27Z Limb of the Sun crests the horizon
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

    from skyfield.api import Star

    galactic_center = Star(ra_hours=(17, 45, 40.04),
                           dec_degrees=(-29, 0, 28.1))

    f = almanac.risings_and_settings(eph, galactic_center, bluffton)
    t, y = almanac.find_discrete(t0, t1, f)

    for ti, yi in zip(t, y):
        verb = 'rises above' if yi else 'sets below'
        print(ti.utc_iso(' '), '- Galactic Center', verb, 'the horizon')

.. testoutput::

    2020-02-01 10:29:00Z - Galactic Center rises above the horizon
    2020-02-01 18:45:46Z - Galactic Center sets below the horizon

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

The phases of the Moon are the same for everyone on Earth,
so you don’t need to specify the longitude and latitude of your location.
Simply ask for the current phase of the Moon.
The return value is an angle
where 0° is New Moon, 90° is First Quarter,
180° is Full Moon, and 270° is Last Quarter:

.. testcode::

    t = ts.utc(2020, 11, 19)
    phase = almanac.moon_phase(eph, t)
    print('Moon phase: {:.1f} degrees'.format(phase.degrees))

.. testoutput::

    Moon phase: 51.3 degrees

Or you can have Skyfield search over a range of dates for the moments
when the Moon reaches First Quarter, Full Moon, Last Quarter, and New Moon:

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
    [1 0]
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

    ['2019-09-02T10:42:26Z', '2020-10-13T23:25:55Z']
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

.. _transits:

Meridian Transits
=================

Every day the Earth’s rotation
swings the sky through nearly 360°,
leaving the celestial poles stationary
while bringing each star and planet in turn
across your *meridian* —
the line of right ascension in the sky above you
that runs from the South Pole to the North Pole through your local zenith.

You can ask Skyfield for the times at which a body
crosses your meridian:

.. testcode::

    t0 = ts.utc(2020, 11, 6)
    t1 = ts.utc(2020, 11, 8)
    t = almanac.find_transits(observer, eph['Mars'], t0, t1)

    print(t.utc_strftime('%Y-%m-%d %H:%M'))

.. testoutput::

    ['2020-11-06 03:32', '2020-11-07 03:28']

Skyfield also has an older mechanism for detecting transits
that isn’t as fast but that also returns the moments of anti-transit,
when a body crosses the line of right ascension that crosses your local nadir:

.. testcode::

    t0 = ts.utc(2020, 11, 6)
    t1 = ts.utc(2020, 11, 7)
    f = almanac.meridian_transits(eph, eph['Mars'], bluffton)
    t, y = almanac.find_discrete(t0, t1, f)

    print(t.utc_strftime('%Y-%m-%d %H:%M'))
    print(y)
    print([almanac.MERIDIAN_TRANSITS[yi] for yi in y])

.. testoutput::

    ['2020-11-06 03:32', '2020-11-06 15:30']
    [1 0]
    ['Meridian transit', 'Antimeridian transit']

Some astronomers call these moments
“upper culmination” and “lower culmination” instead.

Observers often think of transit as the moment
when an object is highest in the sky,
but that’s only roughly true.
At very high precision,
if the body has any north or south velocity
then its moment of highest altitude will be slightly earlier or later.

Bodies near the poles are exceptions to the general rule
that a body is visible at transit but below the horizon at antitransit.
For a body that’s circumpolar from your location,
transit and antitransit are both moments of visibility,
when it stands above and below the pole.
And objects close to the opposite pole will always be below the horizon,
even as they invisibly transit your line of longitude
down below your horizon.

TODO
====

In case you are maintaining older code,
versions of Skyfield before 1.47 could only compute sunrises and sunsets
with an almanac routine
that was both slower than the routine described above,
and that also tended to miss sunrises and sunsets in the Arctic and Antarctic.
Here’s how the older routine is called:

.. testcode::

    t0 = ts.utc(2018, 9, 12, 4)
    t1 = ts.utc(2018, 9, 13, 4)
    t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(eph, bluffton))

    print(t.utc_iso())
    print(y)

.. testoutput::

    ['2018-09-12T11:13:13Z', '2018-09-12T23:49:38Z']
    [1 0]

The result ``t`` will be an array of times, and ``y`` will be ``1`` if
the sun rises at the corresponding time and ``0`` if it sets.

If you need to provide your own custom value for refraction, adjust the
estimate of the Sun’s radius, or account for a vantage point above the
Earth’s surface, see :ref:`risings-and-settings` to learn about the more
versatile :func:`~skyfield.almanac.risings_and_settings()` routine.

Moonrise and moonset
====================

Passing the Moon to the
:func:`~skyfield.almanac.find_risings()`
and
:func:`~skyfield.almanac.find_settings()`
routines that were discussed in the previous section
will let you find the times of moonrise and moonset.
Skyfield uses the
`official definition of moonrise and moonset
<https://aa.usno.navy.mil/faq/RST_defs>`_
from the United States Naval Observatory:
the moment when the top edge of the Moon
is exactly 34 arcminutes below the horizon,
as an approximate correction for atmospheric refraction.

.. testcode::

    moon = eph['Moon']
    t0 = ts.utc(2023, 12, 27)
    t1 = ts.utc(2023, 12, 29)

    t, y = almanac.find_risings(observer, moon, t0, t1)
    print('Moonrises (UTC):', t.utc_iso(' '))

    t, y = almanac.find_settings(observer, moon, t0, t1)
    print('Moonsets (UTC):', t.utc_iso(' '))

.. testoutput::

    Moonrises (UTC): ['2023-12-27 22:40:11Z', '2023-12-28 23:43:48Z']
    Moonsets (UTC): ['2023-12-27 13:54:47Z', '2023-12-28 14:39:33Z']

Read the previous section to learn about the Boolean array ``y``,
which is ``False`` for Arctic and Antarctic locations
when the Moon reaches the meridian without crossing the horizon.

Twilight
========

The routine :func:`~skyfield.almanac.dark_twilight_day()`
returns a separate code for each of the phases of twilight:

0. Dark of night.
1. Astronomical twilight.
2. Nautical twilight.
3. Civil twilight.
4. Daytime.

You can find a full example of its use
at the :ref:`dark_twilight_day() example`.

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

.. _lunar-eclipses:

Lunar eclipses
==============

Skyfield can find the dates of lunar eclipses.

.. testcode::

    from skyfield import eclipselib

    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2020, 1, 1)
    t, y, details = eclipselib.lunar_eclipses(t0, t1, eph)

    for ti, yi in zip(t, y):
        print(ti.utc_strftime('%Y-%m-%d %H:%M'),
              'y={}'.format(yi),
              eclipselib.LUNAR_ECLIPSES[yi])

.. testoutput::

    2019-01-21 05:12 y=2 Total
    2019-07-16 21:31 y=1 Partial

Note that any eclipse forecast
is forced to make arbitrary distinctions
when eclipses fall very close to the boundary
between the categories “partial”, “penumbral”, and “total”.
Skyfield searches for lunar eclipses using the techniques described
in the *Explanatory Supplement to the Astronomical Almanac.*
Here is its current behavior:

.. Note to myself: these claims are generated by editing and re-running
   the ./design/eclipses_lunar.py script.

* Skyfield currently finds every one of the 3,642 lunar eclipses
  listed for the years AD 1000–2500
  in NASA’s
  `Five Millennium Canon of Lunar Eclipses
  <https://eclipse.gsfc.nasa.gov/SEpubs/5MCLE.html>`_
  by Espenak and Meeus.

* But some slight disagreements are inevitable,
  because Skyfield uses a modern ephemeris for Earth and Moon positions,
  while the *Supplement* used the old VSOP87 theory.
  In 8 cases over the years AD 1000–2500 (around 0.2% of the eclipses listed),
  Skyfield disagrees with the *Canon*
  about whether an eclipse was partial or total.
  And on 1571 July 7 Skyfield finds an eclipse,
  but the *Canon* judges the Moon to have narrowly missed our shadow
  on that occasion.

* Skyfield tends to return eclipse times
  that are a few seconds earlier than those given by the *Canon*.
  For decades near the present the disagreement
  rarely exceeds 2 seconds,
  but for eclipses 2,000 years ago the difference
  can be as large as 20 seconds.

* Over the full five millennia covered by the *Canon*,
  Skyfield misses only four eclipses, finds two extra eclipses,
  and agrees with the *Canon*\ ’s category
  (partial, penumbral, total)
  more than 99.8% of the time.
  Of the two missing eclipses that are closest to the modern day,
  the *Canon* gives the April 859 eclipse
  a penumbral magnitude of only 0.0007,
  and the February 2791 eclipse
  a penumbral magnitude of only 0.0006 —
  so the missing eclipses were not exactly major celestial events.

To help you study each eclipse in greater detail,
Skyfield returns a ``details`` dictionary of extra arrays
that provide the dimensions of the Moon and of the Earth’s shadow
at the height of the eclipse.
The means of each field is hopefully self-explanatory;
if any of the terms is unfamiliar,
try looking it up online.

.. testcode::

    for name, values in sorted(details.items()):
        print(f'{name:24}  {values}')

.. testoutput::

    closest_approach_radians  [0.00657921 0.01029097]
    moon_radius_radians       [0.00485608 0.00435481]
    penumbra_radius_radians   [0.02278213 0.02077108]
    penumbral_magnitude       [2.16831186 1.70327942]
    umbra_radius_radians      [0.01332129 0.01161176]
    umbral_magnitude          [1.19418911 0.65164729]

The first element in each of these sequences
corresponds to the first eclipse we discovered above, on 2019-01-21,
while the second element belongs to the eclipse on 2019-07-16.

By combining these dimensions
with the position of the Moon at the height of the eclipse
(which you can generate using Skyfield’s usual approach
to computing a position),
you should be able to produce a detailed diagram of each eclipse.

For a review of the parameters that differ between eclipse forecasts,
see NASA’s
`Enlargement of Earth's shadows
<https://eclipse.gsfc.nasa.gov/LEcat5/shadow.html>`_
page on their Five Millennium Canon site.
If you need lunar eclipse forecasts
generated by a very specific set of parameters,
try cutting and pasting Skyfield’s ``lunar_eclipses()`` function
into your own code
and making your adjustments there —
you will have complete control of the outcome,
and your application will be immune
to any tweaking that takes place in Skyfield in the future
if it’s found that Skyfield’s eclipse accuracy can become even better.
