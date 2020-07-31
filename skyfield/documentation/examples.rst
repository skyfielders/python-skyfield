
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
    from skyfield.api import Topos, load

    # Figure out local midnight.
    zone = timezone('US/Eastern')
    now = zone.localize(dt.datetime.now())
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    next_midnight = midnight + dt.timedelta(days=1)

    ts = load.timescale()
    t0 = ts.from_datetime(midnight)
    t1 = ts.from_datetime(next_midnight)
    eph = load('de421.bsp')
    bluffton = Topos('40.8939 N', '83.8917 W')
    f = almanac.dark_twilight_day(eph, bluffton)
    times, events = almanac.find_discrete(t0, t1, f)

    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')

.. testoutput::

    2020-04-19 05:09   Astronomical twilight starts
    2020-04-19 05:46   Nautical twilight starts
    2020-04-19 06:21   Civil twilight starts
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

    ts = load.timescale()
    t = ts.utc(2019, 12, 9, 15, 36)

    eph = load('de421.bsp')
    sun, moon, earth = eph['sun'], eph['moon'], eph['earth']

    e = earth.at(t)
    _, slon, _ = e.observe(sun).apparent().ecliptic_latlon()
    _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon()
    phase = (mlon.degrees - slon.degrees) % 360.0

    print('{0:.1f}'.format(phase))

.. testoutput::

    149.4

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

    from skyfield.api import load, Topos
    from skyfield.trigonometry import position_angle_of

    ts = load.timescale()
    t = ts.utc(2019, 9, 30, 23)

    eph = load('de421.bsp')
    sun, moon, earth = eph['sun'], eph['moon'], eph['earth']
    boston = earth + Topos('42.3583 N', '71.0636 W')

    b = boston.at(t)
    m = b.observe(moon).apparent()
    s = b.observe(sun).apparent()
    print(position_angle_of(m.altaz(), s.altaz()))

.. testoutput::

    238deg 55' 55.3"

The :func:`~skyfield.trigonometry.position_angle_of()` routine
will not only accept the output of :meth:`~skyfield.positionlib.Apparent.altaz()`,
but also of :meth:`~skyfield.positionlib.ICRF.ecliptic_latlon()`
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

    from skyfield.api import Star, Topos, load
    from skyfield.almanac import find_discrete, risings_and_settings
    from pytz import timezone

    ts = load.timescale()
    t0 = ts.utc(2019, 1, 19)
    t1 = ts.utc(2019, 1, 21)

    moab = Topos('38.5725 N', '109.54972238 W')
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
    topos = api.Topos(latitude_degrees=42, longitude_degrees=-87)
    observer = topos.at(t)
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
Calling the position’s
:meth:`~skyfield.positionlib.Geocentric.subpoint()` method
lets you compute the Earth latitude and longitude
from which the position is directly overhead.

But sometimes the right ascension and declination of the position
are known already.
Instead of creating a :class:`~skyfield.starlib.Star` with those coordinates
and asking it to compute its position,
there is a simpler approach:
creating the position directly.

.. testcode::

    from skyfield.api import load
    from skyfield.positionlib import position_of_radec

    ts = load.timescale()
    t = ts.utc(2020, 1, 3, 12, 45)

    earth = 399  # NAIF code for the Earth center of mass
    ra_hours = 3.79
    dec_degrees = 24.1167
    pleiades = position_of_radec(ra_hours, dec_degrees, t=t, center=earth)
    subpoint = pleiades.subpoint()

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

   from skyfield.api import Topos, load
   from skyfield.functions import length_of

   ts = load.timescale()
   t = ts.utc(2019, 1, 1)

   bierstadt = Topos('39.5828 N', '105.6686 W', elevation_m=4287.012)
   m1 = length_of(bierstadt.at(t).position.m)
   print(int(m1))

   accra = Topos('5.6037 N', '0.1870 W', elevation_m=61)
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
