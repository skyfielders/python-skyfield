
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

At what angle is the sun to the crescent Moon?
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

    238deg 55' 57.0"

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

When is the Galactic Center above the horizon?
==============================================

This approach is not specific to the galactic center,
which is used here merely as an example target.
You can use any Skyfield body
that generates a position when passed as an argument
to a call like::

    observer.at(t).observe(...)

Possible targets include the Sun, Moon, planets,
and any star or other distant object
whose right ascension and declination
you use to build a ``Star()`` object.

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

    Sat 19 05:53 MST rises
    Sat 19 14:26 MST sets
    Sun 20 05:49 MST rises
    Sun 20 14:22 MST sets

Which geographic location is farther from Earth’s center?
=========================================================

After my hike of Mount Bierstadt in Colorado,
a friend suggested that its 14,000 feet of elevation
might have carried me farther from the Earth’s center
than I had ever traveled before.
It was a romantic thought,
that under my own power
I had hiked farther from my home planet’s core
than ever I had stood before.

But there was a problem:
I knew that I had once visited a location
only a few degrees away from the Earth’s equator,
and that the Earth’s equatorial bulge
might push even modest elevations at that latitude
out farther from the Earth’s center
than even a mountaintop in Colorado.

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

