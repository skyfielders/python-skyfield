
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

