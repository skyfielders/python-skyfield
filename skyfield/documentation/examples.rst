
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

Here’s the list of examples that you can find in the sections below:

.. contents::
   :local:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup(
       (2020, 4, 19, 17, 58))

And the following sections are the examples themselves.

What time is solar noon, when the Sun transits the meridian?
============================================================

The Earth travels fastest in its orbit
when closest to the Sun in January,
and slowest in July when they are farthest apart.
This stretches out January days by around 10 seconds,
shortens July days by the same amount,
and means that your watch never reads exactly 12:00
when the Sun transits your meridian at “solar noon”.

You can compute solar noon
by asking at what time the Sun transits the meridian at your location.
To learn the time of solar noon today, for example, you might do this:

.. testcode::

    import datetime as dt
    from pytz import timezone
    from skyfield import almanac
    from skyfield.api import N, W, wgs84, load

    zone = timezone('US/Eastern')
    now = zone.localize(dt.datetime.now())
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    next_midnight = midnight + dt.timedelta(days=1)

    ts = load.timescale()
    t0 = ts.from_datetime(midnight)
    t1 = ts.from_datetime(next_midnight)
    eph = load('de421.bsp')
    bluffton = wgs84.latlon(40.8939 * N, 83.8917 * W)

    f = almanac.meridian_transits(eph, eph['Sun'], bluffton)
    times, events = almanac.find_discrete(t0, t1, f)

    # Select transits instead of antitransits.
    times = times[events == 1]

    t = times[0]
    tstr = str(t.astimezone(zone))[:19]
    print('Solar noon:', tstr)

.. testoutput::

    Solar noon: 2020-04-19 13:34:33

.. _dark_twilight_day() example:

When will it get dark tonight?
==============================

Sunrise, sunset, and the several varieties of twilight
are all available through the :doc:`almanac` module.
Here’s the script I use when I want to know when it will be dark enough
to see the stars —
or how early I need to rise to see the morning sky:

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

    previous_e = f(t0).item()
    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        if previous_e < e:
            print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')
        else:
            print(tstr, ' ', almanac.TWILIGHTS[previous_e], 'ends')
        previous_e = e

.. testoutput::

    2020-04-19 05:09   Astronomical twilight starts
    2020-04-19 05:46   Nautical twilight starts
    2020-04-19 06:20   Civil twilight starts
    2020-04-19 06:49   Day starts
    2020-04-19 20:20   Day ends
    2020-04-19 20:48   Civil twilight ends
    2020-04-19 21:23   Nautical twilight ends
    2020-04-19 22:00   Astronomical twilight ends

As you can see from the above code,
if the new light level is brighter
then we say that the new level “starts”,
but if the new level is darker
then we say the previous level “ends” —
so instead of saying “astronomical twilight *starts* at 21:23”
we say “nautical twilight *ends* at 21:23.”
That’s why the code keeps up with ``previous_e``
and compares it to the new level of twilight.

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
Skyfield also has a method for computing
what fraction of a spherical body is illuminated by the sun.

.. testcode::

    from skyfield.api import load
    from skyfield.framelib import ecliptic_frame

    ts = load.timescale()
    t = ts.utc(2019, 12, 9, 15, 36)

    eph = load('de421.bsp')
    sun, moon, earth = eph['sun'], eph['moon'], eph['earth']

    e = earth.at(t)
    s = e.observe(sun).apparent()
    m = e.observe(moon).apparent()

    _, slon, _ = s.frame_latlon(ecliptic_frame)
    _, mlon, _ = m.frame_latlon(ecliptic_frame)
    phase = (mlon.degrees - slon.degrees) % 360.0

    percent = 100.0 * m.fraction_illuminated(sun)

    print('Phase (0°–360°): {0:.1f}'.format(phase))
    print('Percent illuminated: {0:.1f}%'.format(percent))

.. testoutput::

    Phase (0°–360°): 149.4
    Percent illuminated: 92.9%

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

Are planets separated by 0° at conjunction and 180° at opposition?
==================================================================

It surprises many newcomers to astronomy
that the angular separation between two planets
never drops all the way to 0° at conjunction
nor increases all the way to a full 180° at opposition.
The reason is that the planets will still have at least a slight
lingering difference in latitude,
even when their longitudes have brought them together
or have placed them opposite each other in the sky.

We can take as an example
the date and time of the conjunction between Mars and the Sun
computed in the :ref:`oppositions-conjunctions` section of the Almanac page.
How close are they in the sky at that moment?
The :func:`~skyfield.positionlib.ICRF.separation_from()` method
measures raw angular distance
between any two points in the sky:

.. testcode::

    from skyfield.api import load
    from skyfield.framelib import ecliptic_frame

    ts = load.timescale()
    eph = load('de421.bsp')
    sun, mars = eph['sun'], eph['mars']

    t = ts.utc(2019, 9, 2, 10, 42, 26)
    e = earth.at(t)
    s = e.observe(sun).apparent()
    m = e.observe(mars).apparent()
    print('{:.5f}°'.format(m.separation_from(s).degrees))

.. testoutput::

    1.08256°

They are more than one degree apart!
How can that be,
if their ecliptic longitudes are at that moment the same?
Let’s use Skyfield’s :data:`~skyfield.framelib.ecliptic_frame`
to express their positions in :ref:`ecliptic-coordinates`:

.. testcode::

    print('     Latitude Longitude')

    lat, lon, distance = s.frame_latlon(ecliptic_frame)
    print('Sun  {:.5f}° {:.5f}°'.format(lat.degrees, lon.degrees))

    lat, lon, distance = m.frame_latlon(ecliptic_frame)
    print('Mars {:.5f}° {:.5f}°'.format(lat.degrees, lon.degrees))

.. testoutput::

         Latitude Longitude
    Sun  0.00005° 159.68641°
    Mars 1.08260° 159.68641°

While the Sun sits very close to the ecliptic —
as we would expect, since the ecliptic is defined
as the course the Sun takes around the sky each year —
the inclination of the orbit of Mars has carried it
more than one degree above the ecliptic.
That’s why the :func:`~skyfield.positionlib.ICRF.separation_from()` method
still measured an angle of more than one degree between them.

A similar situation pertains at opposition:

.. testcode::

    t = ts.utc(2020, 10, 13, 23, 25, 55)

    e = earth.at(t)
    s = e.observe(sun).apparent()
    m = e.observe(mars).apparent()

    print('Separation: {:.5f}°'.format(m.separation_from(s).degrees))

    print('')
    print('     Latitude Longitude')

    lat, lon, distance = s.frame_latlon(ecliptic_frame)
    print('Sun  {:.5f}° {:.5f}°'.format(lat.degrees, lon.degrees))

    lat, lon, distance = m.frame_latlon(ecliptic_frame)
    print('Mars {:.5f}° {:.5f}°'.format(lat.degrees, lon.degrees))

.. testoutput::

    Separation: 177.00424°

         Latitude Longitude
    Sun  0.00007° 201.07794°
    Mars -2.99582° 21.07794°

Even though their ecliptic longitudes are 180° apart,
the fact that neither the Sun nor Mars is lying exactly on the ecliptic
means that the :func:`~skyfield.positionlib.ICRF.separation_from()` method
finds that they are not quite 180° apart.

In case you run across the term ‘elongation’
in discussions of conjunctions and oppositions,
it’s shorthand for ‘the angle between a planet and the Sun’ —
and so each of the angular separations printed above can,
more specifically,
be labeled as the ‘elongation of Mars’ on those dates.

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
the output of :meth:`~skyfield.positionlib.ICRF.altaz()`,
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

At what rate is a target moving across the sky?
===============================================

If you are interested in the rate at which a target is moving across the sky,
you can call Skyfield’s
:meth:`~skyfield.positionlib.ICRF.frame_latlon_and_rates()` method
and pass it the frame of reference in which you want the angles measured.
First, compute the target’s position relative to your geographic location:

.. testcode::

    from skyfield.api import load, wgs84, N,S,E,W

    ts = load.timescale()
    t = ts.utc(2021, 2, 3, 0, 0)
    planets = load('de421.bsp')
    earth, mars = planets['earth'], planets['mars']
    topos = wgs84.latlon(35.1844 * N, 111.6535 * W, elevation_m=2099.5)

    a = (earth + topos).at(t).observe(mars).apparent()

In Skyfield, a topocentric location object like ``topos``
is also a reference frame oriented to the
:ref:`location’s horizon and zenith <horizontal-coordinates>`.
So if you pass it to the
:meth:`~skyfield.positionlib.ICRF.frame_latlon_and_rates()` method,
Skyfield will compute the rates at which the altitude and azimuth are changing
as the target moves across the sky:

.. testcode::

    (alt, az, distance,
     alt_rate, az_rate, range_rate) = a.frame_latlon_and_rates(topos)

    print('Alt: {:+.1f} asec/min'.format(alt_rate.arcseconds.per_minute))
    print('Az:  {:+.1f} asec/min'.format(az_rate.arcseconds.per_minute))

.. testoutput::

    Alt: +548.7 asec/min
    Az:  +1586.4 asec/min

.. Why is this a comment? Because it is false. Using equinox-of-date
   combines two motions, that of the body across the heavens, and that
   of the Earth’s pole. Should I re-do this using Hour Angle?

   Or maybe I should do this using a frame that is ICRS / J2000.
   Except that, does Skyfield have such a frame?

   Anyway:

   Or, if you instead want to know how fast the target is moving
   against the background of stars,
   you can pass Skyfield’s built-in
   :data:`~skyfield.framelib.true_equator_and_equinox_of_date` reference frame
   to compute rates of moment in right ascension and declination:

   .. testcode::

       from skyfield import framelib

       teqeq = framelib.true_equator_and_equinox_of_date
       (dec, ra, distance,
        dec_rate, ra_rate, range_rate) = a.frame_latlon_and_rates(teqeq)

       print(f'RA:  {ra_rate.arcseconds.per_hour:+.1f} asec/hr')
       print(f'Dec: {dec_rate.arcseconds.per_hour:+.1f} asec/hr')

   .. testoutput::

       RA:  +78.7 asec/hr
       Dec: +25.6 asec/hr

   Note that, contrary to Skyfield’s usual custom,
   this technique returns declination as the first return value
   instead of returning right ascension first.
   That’s because the
   :meth:`~skyfield.positionlib.ICRF.frame_latlon_and_rates()` method
   always returns the latitude-like coordinate first,
   which measures the target’s angle ±90° above or below the plane,
   and then the longitude-like coordinate second,
   which measures the target’s position 0°–360° around.
   The ``latlon`` in its name can help you remember this.

You can choose other units besides ``arcseconds`` and ``per_minute``.
For the possible numerators see
:class:`~skyfield.units.AngleRate`,
and for the possible denominators see
:class:`~skyfield.units.Rate`.

Computing the range rate-of-change
----------------------------------

Both of the calls above return a ``range_rate``
that is positive if the body is moving away
and negative if the target is moving closer:

.. testcode::

    print('Range rate: {:+.1f} km/s'.format(range_rate.km_per_s))

.. testoutput::

    Range rate: +16.8 km/s

Computing angular speed
-----------------------

You might think that you could compute
a target’s total angular speed across the sky
by simply subjecting the two angular rates of change
to the Pythagorean theorem.

But that won’t work, because of a subtlety:
it turns out that all of the different kinds of longitude —
including right ascension, azimuth, and ecliptic longitude —
have lines that are far apart at the equator
but that draw closer and closer together near the poles.
I hope that there is an elegant antique globe sitting near you
as you read this in your armchair.
Look at the lines on its surface.
Down at the equator,
the lines of longitude stand far apart,
and to move 15° in longitude
you would have to travel across very nearly 15° of the Earth’s surface.
But now look at the poles.
The lines of longitude draw so close together that,
if you’re close enough to the pole,
you could cross 15° of longitude
by traveling only a very short distance!

Happily, spherical trigonometry gives us a simple correction to apply.
Multiplying the longitude rate by the cosine of the latitude
gives a bare angular rate of motion across the sky,
that can safely be tossed into the Pythagorean theorem:

.. testcode::

    from numpy import cos, sqrt

    ralt = alt_rate.degrees.per_minute
    raz = az_rate.degrees.per_minute * cos(alt.radians)

    degrees_per_minute = sqrt(ralt*ralt + raz*raz)
    print('{:.4f}° per minute'.format(degrees_per_minute))

.. testoutput::

    0.2392° per minute

In exactly the same way,
if instead you wanted to compute a target’s speed
against the background of stars,
you would multiply the rate at which the right ascension is changing
by the cosine of the declination
before combining them with the Pythagorean theorem.

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
