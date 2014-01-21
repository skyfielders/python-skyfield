
===========================
 Positions and Coordinates
===========================

.. currentmodule:: skyfield.api

Skyfield is careful to distinguish the *position* of an object
from the several choices of *coordinate*
that can be used to designate that position with numbers.
There are only three basic kinds of position that Skyfield recognizes,
but several different ways in which a position
can be turned into coordinates.

Here is a quick reference to the three basic kinds of position,
together with all of the attributes and methods that they support:

.. parsed-literal::

    Three positions

    obj(time)           →  Barycentric position (BCRS)
     └─ observe(obj2)   →  Astrometric position (ΔBCRS)
         └─ apparent()  →  Apparent position (GCRS)

    Barycentric, Astrometric, or Apparent position
     │
     ├── `position <api.html#Position.position>`_.AU         →   x, y, z
     ├── `position <api.html#Position.position>`_.km         →   x, y, z
     ├── `position.to(unit) <api.html#Distance.to>`_   →   x, y, z
     │
     ├── `velocity <api.html#Position.velocity>`_.AU_per_d   →   xdot, ydot, zdot
     ├── `velocity <api.html#Position.velocity>`_.km_per_s   →   xdot, ydot, zdot
     ├── `velocity.to(unit) <api.html#Distance.to>`_   →   xdot, ydot, zdot
     │
     ├── `radec() <api.html#Position.radec>`_             →   ra, dec, distance
     └── `radec(epoch=jd) <api.html#Position.radec>`_     →   ra, dec, distance

    Apparent position only
     │
     └── `altaz() <api.html#Position.altaz>`_             →   alt, az, distance

    Angle like ra, dec, alt, and az
     │
     ├── `radians() <api.html#Angle.radians>`_           →   6.266029488577352
     │
     ├── `hours() <api.html#Angle.hours>`_             →   23.934469599999996
     ├── `hms() <api.html#Angle.hms>`_               →   (1, 23, 56, 4, 0)
     ├── `hstr() <api.html#Angle.hstr>`_              →   '23h 56m 04.09s'
     ├── `hstr(places=4) <api.html#Angle.hstr>`_      →   '23h 56m 04.0906s'
     │
     ├── `degrees() <api.html#Angle.degrees>`_           →   359.017044
     ├── `dms() <api.html#Angle.dms>`_               →   (1, 359, 1, 1, 0)
     ├── `dstr() <api.html#Angle.dstr>`_              →   '359deg 01\' 01.4"'
     └── `dstr(places=3) <api.html#Angle.dstr>`_      →   '359deg 01\' 01.358"'

The rest of this page is simply designed to explain
all of the features outlined in the quick reference above.
All hyperlinked attributes and method names,
both in the text above and in the explanations below,
lead to the low-level :doc:`api`
which explains each option in even greater detail.

Generating positions
====================

Skyfield already supports two kinds of object
that can compute their position,
and will soon be supporting more.
Each object is a simple Python callable
that can take a date as its argument
and return its position as of that moment:

* :doc:`planets` — the eight planets and Pluto are all supported,
  thanks to the excellent work of the Jet Propulsion Laboratory (JPL)
  and to Skyfield’s support for their major solar system ephemerides. ::

    from skyfield.api import earth, now
    position = earth(now())

* :doc:`earth-satellites` — Earth satellite positions can be generated
  from public TLE elements describing their current orbit,
  which you can download from Celestrak. ::

    from skyfield.api import earth, now
    sat = earth.satellite()
    position = sat(now())

Read :doc:`time` for more information
about how to build dates and pass them to planets and satellites
when generating positions.

What can you do with a position once it has been generated?
The rest of this document is a complete tour of the possibilities.

Barycentric position
====================

When you ask Skyfield for the position of a planet or star,
it produces a three-dimensional position
that is measured from the Solar System *barycenter* —
the center of mass around which all of the planets revolve.
The position is stored as *x*, *y*, and *z* coordinates
in the International Celestial Reference System (ICRS),
a permanent frame of reference
that is a high-precision replacement
for the old J2000.0 system
that was popular at the end of the 20th century.

An ICRS coordinate centered on the solar system barycenter
is called a *Barycentric Celestial Reference System* (BCRS) coordinate.

You can view these coordinates
by asking Skyfield for their :attr:`~Position.position` attribute:

.. testcode::

    # BCRS positions of Earth and Jupiter

    from skyfield.api import earth, jupiter
    print earth(utc=(1980, 1, 1)).position.AU
    print jupiter(utc=(1980, 1, 1)).position.AU

.. testoutput::

    [-0.16287311  0.88787399  0.38473904]
    [-4.71061475  2.32932129  1.11328106]

The coordinates shown above are measured
using the Astronomical Unit (AU),
which is the average distance from the Earth to the Sun.
So the value ``-4.71`` indicates a distance
nearly five times farther from the Sun than that of the Earth.
You can, if you want, ask for these coordinates
in kilometers with the :attr:`~Position.km` attribute.
And if you have the third-party AstroPy package installed,
then you can convert these coordinates
into any length unit with the :meth:`~Position.to()` method.

Astrometric position
====================

You might think that you could determine
the position of Jupiter in the night sky
by simply subtracting these two positions
to generate the vector difference between them.
But that would ignore the fact that light takes several minutes
to travel between Jupiter and the Earth.
The image of Jupiter in our sky
does not show us where it *is*, right now,
but where it *was* — several minutes ago —
when the light now reaching our eyes or instruments
actually left its surface.

Correcting for the light-travel time
does not simply fix a minor inconvenience,
but reflects a very deep physical reality.
Not only the light from Jupiter,
but *all* of its physical effects,
arrive no faster than the speed of light.
As Jupiter tugs us with its gravity,
we do not get pulled in the direction of the “real” Jupiter —
we get tugged in the direction of its time-delayed image
hanging in the sky above us!

So Skyfield offers a :meth:`~Position.observe()` method
that carefully backdates the position of another object
to determine where it was when it generated the image
that we see in our sky:

.. testcode::

    # Observing Jupiter from the Earth's position

    astro = earth(utc=(1980, 1, 1)).observe(jupiter)
    print astro.position.AU

.. testoutput::

    [-4.54763822  1.44160883  0.72860876]

This light-delayed position is called the *astrometric* position,
and is traditionally mapped on a star chart
by the angles *right ascension* and *declination*
that you can compute using the :meth:`~Position.radec()` method
and display using their :meth:`~Angle.hstr()`
and :meth:`~Angle.dstr()` methods:

.. testcode::

    # Astrometric RA and declination

    ra, dec, distance = astro.radec()
    print ra.hstr()
    print dec.dstr()
    print distance.AU

.. testoutput::

    10h 49m 38.71s
    +08deg 41' 00.6"
    4.82598384993

As we will explore in the next section,
objects never appear at exactly the position in the sky
predicted by the simple and ideal astrometric position.
But it is useful for mapping the planet
against the background of stars in a
`printed star atlas <http://www.amazon.com/s/?_encoding=UTF8&camp=1789&creative=390957&linkCode=ur2&pageMinusResults=1&suo=1389754954253&tag=letsdisthemat-20&url=search-alias%3Daps#/ref=nb_sb_noss_1?url=search-alias%3Daps&field-keywords=star%20atlas&sprefix=star+%2Caps&rh=i%3Aaps%2Ck%3Astar%20atlas&sepatfbtf=true&tc=1389754955568>`_,
because star atlases also use astrometric positions.

Apparent position
=================

To determine the position of an object in the night sky
with even greater accuracy,
two further effects must be taken into account:

*Deflection*
  The object’s light is bent,
  and thus its image displaced,
  if the light passes close to another large mass
  on its way to the observer.
  This will happen if the object lies very near to the Sun in the sky,
  for example, or is nearly behind Jupiter.
  The effect is small,
  but must be taken into account for research-grade results.

*Aberration*
  The velocity of the Earth itself through space
  adds a very slight slant to light arriving at our planet,
  in the same way that rain or snow
  seen through the windshield while driving
  appears to be slanting towards you
  because of your own motion.
  The effect is small — at most about 20 arcseconds —
  and so was not discovered until 1729.
  The discovery finally proved that the Earth goes around the Sun.

Skyfield lets you apply both of these effects
by invoking the :meth:`~Astrometric.apparent()` method.
Like an astrometric position, an apparent position
is typically expressed as the angles
*right ascension* and *declination*:

.. testcode::

    # Apparent GCRS ("J2000.0") coordinates

    apparent = astro.apparent()
    ra, dec, distance = apparent.radec()

    print ra.hstr()
    print dec.dstr()
    print distance.AU

.. testoutput::

    10h 49m 39.34s
    +08deg 40' 56.4"
    4.82598384993

But it is actually unusual to print apparent coordinates
in a permanent unchanging reference frame like the ICRS,
so you are unlikely to find the two values above
if you look up the position of Jupiter on 1980 January 1
in an almanac or by using other astronomy software.

Instead, apparent positions are usually expressed
relative to the Earth’s real orientation
as its rolls and tumbles through space —
which, after all,
is how right ascension and declination were defined
through most of human history,
before the invention of the ICRS axes.
The Earth’s poles and equator move at least slightly every day,
and move by very large amounts as years add up to centuries.

To ask for right ascension and declination
relative to the real pole and equator of Earth,
and not the ideal permanent axes of the ICRS,
simply add the keyword argument ``epoch='date'``
when you ask the apparent position for coordinates:

.. testcode::

    # Coordinates relative to true equator and equinox

    ra, dec, distance = apparent.radec(epoch='date')

    print ra.hstr()
    print dec.dstr()
    print distance.AU

.. testoutput::

    10h 48m 36.02s
    +08deg 47' 18.6"
    4.82598384993

These are the coordinates
that should match other astronomy software
and the data in the
`Astronomical Almanac <http://www.amazon.com/s/?_encoding=UTF8&camp=1789&creative=390957&field-keywords=astronomical%20almanac&linkCode=ur2&tag=letsdisthemat-20&url=search-alias%3Daps>`_,
and are sometimes said to be expressed
in the “dynamical reference system” defined by the Earth itself.

Azimuth and altitude
====================

The final result that many users seek
is the *altitude* and *azimuth* of an object
relative to their own location on the Earth’s surface.

The altitude measures the angle above or below the horizon,
with a positive number of degrees meaning “above”
and a negative number indicating that the object
is below the horizon (and impossible to view).
Azimuth measures the angle around the sky from the north pole,
so 0° means that the object is straight north,
90° indicates that the object lies to the east,
180° means south, and 270° means that the object is straight west.

Altitude and azimuth are computed
by calling the :meth:`~Apparent.altaz()` method on an apparent position.
But because the method needs to know whose local horizon to use,
it does not work
on the plain geocentric (“Earth centered”) positions
that have been generating so far:

.. testcode::

    alt, az, distance = apparent.altaz()

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: to compute an apparent position, you must observe from a specific Earth location that you specify using a Topos instance

Instead, you have to give Skyfield your geographic location.
Astronomers use the term *topocentric*
for a position measured relative to a specific location on Earth,
so Skyfield represents Earth locations using a :class:`Topos` class
that you can generate by using the :meth:`Earth.topos` method
of an Earth object:

.. testcode::

    # Altitude and azimuth in the sky of a
    # specific geographic location

    boston = earth.topos('42.3583 N', '71.0603 W')
    astro = boston(utc=(1980, 1, 1)).observe(jupiter)
    alt, az, distance = astro.apparent().altaz()
    print alt.dstr()
    print az.dstr()
    print distance.AU

.. testoutput::

    -23deg 22' 47.8"
    51deg 43' 29.6"
    4.82600082194

So Jupiter was more than 23° below the horizon for Bostonians
on 1980 January 1 at midnight UTC.
