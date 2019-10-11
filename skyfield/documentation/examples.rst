
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

When is the Galactic Center above the horizon?
==============================================

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

