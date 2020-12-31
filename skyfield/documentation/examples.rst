
==========
 Examples
==========

Many programmers work most efficiently
when they can start with a working example
and adapt it to their needs.
Each of the following examples
pulls together several Skyfield features
to solve a general problem,
that should provide readers with a basis
for solving other similar problems of their own.

When will it get dark tonight?
==============================

Sunrise, sunset, and the several varieties of twilight
are all available through the :doc:`almanac` module.
Here’s the script I use when I want to know when it will be dark enough
to see the stars —
or how early I need to rise to see the morning sky:

.. TODO Figure out how to use the timezone itself to find day start/end.

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup(
       (2020, 4, 19, 17, 58))

.. testcode::

    import datetime as dt
    from pytz import timezone
    from skyfield import almanac
    from skyfield.api import N, W, wgs84, load

    # Figure out local midnight.
    zone = timezone('US/Eastern')
    now = zone.localize(dt.datetime.now())
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    next_midnight = midnight + dt.timedelta(days=1)

    ts = load.timescale()
    t0 = ts.from_datetime(midnight)
    t1 = ts.from_datetime(next_midnight)
    eph = load('de421.bsp')
    bluffton = wgs84.latlon(40.8939 * N, 83.8917 * W)
    f = almanac.dark_twilight_day(eph, bluffton)
    times, events = almanac.find_discrete(t0, t1, f)

    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')

.. testoutput::

    2020-04-19 05:09   Astronomical twilight starts
    2020-04-19 05:46   Nautical twilight starts
    2020-04-19 06:20   Civil twilight starts
    2020-04-19 06:49   Day starts
    2020-04-19 20:20   Civil twilight starts
    2020-04-19 20:48   Nautical twilight starts
    2020-04-19 21:23   Astronomical twilight starts
    2020-04-19 22:00   Night starts

What phase is the Moon tonight?
===============================

The *phase* of the Moon is defined
as the angle between the Moon and the Sun along the ecliptic.
This angle is computed as the difference in the *ecliptic longitude*
of the Moon and of the Sun.
The result is an angle that is 0° for the New Moon,
90° at the First Quarter,
180° at the Full Moon,
and 270° at the Last Quarter.

.. testcode::

    from skyfield.api import load
    from skyfield.framelib import ecliptic_frame

    ts = load.timescale()
    t = ts.utc(2019, 12, 9, 15, 36)

    eph = load('de421.bsp')
    sun, moon, earth = eph['sun'], eph['moon'], eph['earth']

    e = earth.at(t)
    _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
    _, mlon, _ = e.observe(moon).apparent().frame_latlon(ecliptic_frame)
    phase = (mlon.degrees - slon.degrees) % 360.0

    print('{0:.1f}'.format(phase))

.. testoutput::

    149.4

What is the angular diameter of a planet, given its radius?
===========================================================

Be careful to select the correct radius
when predicting a planet’s angular diameter in the sky.
Many web sites will quote some kind of “mean radius”
that averages between a planet’s squat polar radius
and its wide equatorial radius.
But most astronomers instead want to know the maximum, not average, diameter
across a planet’s visible face —
so you will want to use the planet’s equatorial radius in your calculation.

For example, a good current estimate of Neptune’s equatorial radius
is 24,764 km.
We would therefore predicts its angular diameter as:

.. testcode::

    import numpy as np
    from skyfield.api import Angle, load

    ts = load.timescale()
    time = ts.utc(2020, 12, 30)

    eph = load('de421.bsp')
    earth, neptune = eph['earth'], eph['neptune barycenter']
    radius_km = 24764.0

    astrometric = earth.at(time).observe(neptune)
    ra, dec, distance = astrometric.apparent().radec()
    apparent_diameter = Angle(radians=np.arcsin(radius_km / distance.km) * 2.0)
    print('{:.6f} arcseconds'.format(apparent_diameter.arcseconds()))

.. testoutput::

    2.257190 arcseconds

This agrees exactly with the output of the NASA HORIZONS system.

When is Venus at its greatest east and west elongations from the Sun?
=====================================================================

This example illustrates the several practical steps
that are often required to both find events of interest
and then to learn more details about them.

* The concept of “elongation from the Sun” is here explained to Skyfield
  with a function that for any given time ``t``
  returns the separation between the Sun and Venus in the sky.

* The ``find_maxima()`` routine is then set loose to find the moments
  over the 3 years 2019–2021 at which the elongation of Venus from the Sun
  reaches its maximum values.
  Skyfield starts by computing the elongation every ``step_days`` = 15 days
  between the search’s start time and end time,
  then hones in everywhere it sees a local maximum:
  a value that’s bigger than either of the two values next to it.

* Finally, a ``for`` loop over the results not only displays each maximum
  but computes and displays an extra fact:
  whether the elongation is an east or west maximum elongation,
  which is defined as whether Venus’s apparent ecliptic longitude
  is ahead of or behind the Sun’s along the great circle of the ecliptic.

This example can serve as a template for many other kinds of custom search:

.. testcode::

    from skyfield.api import load
    from skyfield.framelib import ecliptic_frame
    from skyfield.searchlib import find_maxima

    ts = load.timescale()
    t0 = ts.utc(2019)
    t1 = ts.utc(2022)

    eph = load('de421.bsp')
    sun, earth, venus = eph['sun'], eph['earth'], eph['venus']

    def elongation_at(t):
        e = earth.at(t)
        s = e.observe(sun).apparent()
        v = e.observe(venus).apparent()
        return s.separation_from(v).degrees

    elongation_at.step_days = 15.0

    times, elongations = find_maxima(t0, t1, elongation_at)

    for t, elongation_degrees in zip(times, elongations):
        e = earth.at(t)
        _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
        _, vlon, _ = e.observe(venus).apparent().frame_latlon(ecliptic_frame)
        is_east = (vlon.degrees - slon.degrees) % 360.0 < 180.0
        direction = 'east' if is_east else 'west'
        print('{}  {:4.1f}° {} elongation'.format(
            t.utc_strftime(), elongation_degrees, direction))

.. testoutput::

    2019-01-06 04:53:35 UTC  47.0° west elongation
    2020-03-24 22:13:32 UTC  46.1° east elongation
    2020-08-13 00:14:12 UTC  45.8° west elongation
    2021-10-29 20:51:56 UTC  47.0° east elongation

At what angle in the sky is the crescent Moon?
==============================================

The angle of the crescent Moon changes with the seasons.
In the spring,
a crescent Moon will stand high above the Sun
and appear to be lit from below.
In the autumn,
the Moon sets farther from the Sun along the horizon
and is illuminated more from the side.
What if we wanted to know the exact angle?

You can find the answer
by asking for the Sun’s “position angle” relative to the Moon,
an angle you can compute between any two Skyfield positions.
The angle will be 90° if the Sun is left of the moon,
180° if the Sun is directly below,
and 270° if the Sun is to the right of the Moon.

.. testcode::

    from skyfield.api import N, W, load, wgs84
    from skyfield.trigonometry import position_angle_of

    ts = load.timescale()
    t = ts.utc(2019, 9, 30, 23)

    eph = load('de421.bsp')
    sun, moon, earth = eph['sun'], eph['moon'], eph['earth']
    boston = earth + wgs84.latlon(42.3583 * N, 71.0636 * W)

    b = boston.at(t)
    m = b.observe(moon).apparent()
    s = b.observe(sun).apparent()
    print(position_angle_of(m.altaz(), s.altaz()))

.. testoutput::

    238deg 55' 55.3"

The :func:`~skyfield.trigonometry.position_angle_of()` routine
will not only accept
the output of :meth:`~skyfield.positionlib.Apparent.altaz()`,
but also of :meth:`~skyfield.positionlib.ICRF.frame_latlon()`
if you want a position angle relative to the ecliptic’s north pole.

Beware, though, that :meth:`~skyfield.positionlib.ICRF.radec()`
produces coordinates in the opposite order
from what :func:`~skyfield.trigonometry.position_angle_of()` expects:
right ascension is like longitude, not latitude.
Try reversing the coordinates, like:

.. testcode::

    print(position_angle_of(m.radec(), s.radec()))

.. testoutput::

    282deg 28' 15.7"

Drat, but this angle is backwards, because right ascension increases
toward the east whereas the other angles, like azimuth, increase the
other way around the circle.

When is a body or fixed coordinate above the horizon?
=====================================================

The following code will determine
when the Galactic Center is above the horizon.
The Galactic Center is an example of a fixed object,
like a star or nebula or galaxy,
whose right ascension and declination can be plugged in to a ``Star()`` object.
The code will also work with a body from an ephemeris,
like the Sun, Moon, or one of the planets.

.. testcode::

    from skyfield.api import N, Star, W, wgs84, load
    from skyfield.almanac import find_discrete, risings_and_settings
    from pytz import timezone

    ts = load.timescale()
    t0 = ts.utc(2019, 1, 19)
    t1 = ts.utc(2019, 1, 21)

    moab = wgs84.latlon(38.5725 * N, 109.54972238 * W)
    eph = load('de421.bsp')
    gc = Star(ra_hours=(17, 45, 40.04), dec_degrees=(-29, 0, 28.1))

    f = risings_and_settings(eph, gc, moab)
    tz = timezone('US/Mountain')

    for t, updown in zip(*find_discrete(t0, t1, f)):
        print(t.astimezone(tz).strftime('%a %d %H:%M'), 'MST',
              'rises' if updown else 'sets')

.. testoutput::

    Sat 19 05:51 MST rises
    Sat 19 14:27 MST sets
    Sun 20 05:47 MST rises
    Sun 20 14:23 MST sets

What is the right ascension and declination of a point in the sky?
==================================================================

An observer is often interested in the astronomical coordinates
of a particular position in the sky above them.
If the observer can specify the position
using altitude and azimuth coordinates,
then Skyfield can return its right ascension and declination.

.. testcode::

    from skyfield import api

    ts = api.load.timescale()
    t = ts.utc(2019, 9, 13, 20)
    geographic = api.wgs84.latlon(latitude_degrees=42, longitude_degrees=-87)
    observer = geographic.at(t)
    pos = observer.from_altaz(alt_degrees=90, az_degrees=0)

    ra, dec, distance = pos.radec()
    print(ra)
    print(dec)

.. testoutput::

    13h 41m 14.65s
    +42deg 05' 50.0"

What latitude and longitude is beneath this right ascension and declination?
============================================================================

Most Skyfield calculations,
like an observation of a planet or an Earth satellite,
directly produce a vector position centered on the Earth.
You can pass such a vector
to the :meth:`~skyfield.toposlib.Geoid.subpoint()` method
of a standard geoid to compute latitude and longitude.

But sometimes the right ascension and declination of the position
are known already.
Instead of creating a :class:`~skyfield.starlib.Star` with those coordinates
and asking it to compute its position,
there is a simpler approach:
creating the position directly.

.. testcode::

    from skyfield.api import load, wgs84
    from skyfield.positionlib import position_of_radec

    ts = load.timescale()
    t = ts.utc(2020, 1, 3, 12, 45)

    earth = 399  # NAIF code for the Earth center of mass
    ra_hours = 3.79
    dec_degrees = 24.1167
    pleiades = position_of_radec(ra_hours, dec_degrees, t=t, center=earth)
    subpoint = wgs84.subpoint(pleiades)

    print('Latitude:', subpoint.latitude)
    print('Longitude:', subpoint.longitude)

.. testoutput::

    Latitude: 24deg 10' 33.5"
    Longitude: 123deg 16' 53.9"

Which geographic location is farther from Earth’s center?
=========================================================

After I hiked Mount Bierstadt in Colorado,
a friend suggested that its 14,000 feet of elevation
might have carried me farther from the Earth’s center
than I had ever traveled before.
It was a romantic thought:
that under my own power
I had hiked farther from my home planet’s core
than ever before.

But there was a problem.
I knew that I had once visited a city
only a few degrees away from the Earth’s equator,
and that the Earth’s equatorial bulge
might push even modest elevations at that latitude
out farther from the Earth’s center
than a mountaintop in Colorado.

So I wrote a quick Skyfield script
to compare the distance from the Earth’s center
to both Accra, Ghana, and the top of Mount Bierstadt in Colorado.

.. testcode::

   from skyfield.api import N, W, wgs84, load
   from skyfield.functions import length_of

   ts = load.timescale()
   t = ts.utc(2019, 1, 1)

   bierstadt = wgs84.latlon(39.5828 * N, 105.6686 * W, elevation_m=4287.012)
   m1 = length_of(bierstadt.at(t).position.m)
   print(int(m1))

   accra = wgs84.latlon(5.6037 * N, 0.1870 * W, elevation_m=61)
   m2 = length_of(accra.at(t).position.m)
   print(int(m2))

   assert m2 > m1
   print("I was", int(m2 - m1), "meters farther from the Earth's center\n"
         "when I visited Accra, at nearly sea level, than atop\n"
         "Mt. Bierstadt in Colorado.")

.. testoutput::

    6373784
    6377995
    I was 4211 meters farther from the Earth's center
    when I visited Accra, at nearly sea level, than atop
    Mt. Bierstadt in Colorado.

.. testcleanup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()
