
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
     ├── `position <api.html#Position.position>`_.au         →   x, y, z
     ├── `position <api.html#Position.position>`_.km         →   x, y, z
     ├── `position.to(unit) <api.html#Distance.to>`_   →   x, y, z
     │
     ├── `velocity <api.html#Position.velocity>`_.au_per_d   →   xdot, ydot, zdot
     ├── `velocity <api.html#Position.velocity>`_.km_per_s   →   xdot, ydot, zdot
     ├── `velocity.to(unit) <api.html#Distance.to>`_   →   xdot, ydot, zdot
     │
     ├── `radec(epoch=jd) <api.html#Position.radec>`_     →   ra, dec, distance
     ├── `radec() <api.html#Position.radec>`_             →   ra, dec, distance
     └── `distance() <api.html#Position.distance>`_          →   distance

    Apparent position only
     │
     └── `altaz(…) <api.html#Position.altaz>`_            →   alt, az, distance

    Angle like ra, dec, alt, and az
     │
     ├── `radians <api.html#Angle.radians>`_             →   6.266029488577352
     │
     ├── `hours <api.html#Angle.hours>`_               →   23.934469599999996
     ├── `hstr() <api.html#Angle.hstr>`_              →   '23h 56m 04.09s'
     ├── `hstr(places=4) <api.html#Angle.hstr>`_      →   '23h 56m 04.0906s'
     ├── `hms() <api.html#Angle.hms>`_               →   (23.0, 56.0, 4.0)
     ├── `signed_hms() <api.html#Angle.hms>`_        →   (1.0, 23.0, 56.0, 4.0)
     │
     ├── `degrees <api.html#Angle.degrees>`_             →   359.017044
     ├── `dstr() <api.html#Angle.dstr>`_              →   '359deg 01\' 01.4"'
     ├── `dstr(places=3) <api.html#Angle.dstr>`_      →   '359deg 01\' 01.358"'
     └── `signed_dms() <api.html#Angle.dms>`_        →   (1.0, 359.0, 1.0, 1.0)

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
that can take either a :doc:`Julian date <time>`
or a :ref:`whole Julian date array <date-arrays>`
as its argument and return a corresponding number of positions.

**The planets**
  The eight planets and Pluto are all supported,
  thanks to the excellent work of the Jet Propulsion Laboratory (JPL)
  and Skyfield’s support for their major solar system ephemerides.
  :doc:`Read more <planets>`

  .. testcode::

    from skyfield.api import load, now

    jd = now()
    planets = load('de421.bsp')
    boston = planets('earth').topos('42.3583 N', '71.0603 W')

    # Geocentric

    barycentric = planets('mars').at(jd)
    astrometric = planets('earth').observe('mars').at(jd)
    apparent = planets('earth').observe('mars').at(jd).apparent()

    # Topocentric

    astrometric = boston.observe('mars').at(jd)
    apparent = boston.observe('mars').at(jd).apparent()

**The stars**
  Stars and other fixed objects with catalog coordinates
  generate their current astrometric position
  when observed from a planet. :doc:`Read more <stars>`

  .. TODO - turn the following back into test code

  ::

    from skyfield.api import Star, now

    jd = now()
    boston = earth.topos('42.3583 N', '71.0603 W')
    barnard = Star(ra_hours=(17, 57, 48.49803),
                   dec_degrees=(4, 41, 36.2072))

    # Geocentric

    astrometric = earth(jd).observe(barnard)
    apparent = earth(jd).observe(barnard).apparent()

    # Topocentric

    astrometric = boston(jd).observe(barnard)
    apparent = boston(jd).observe(barnard).apparent()

**Earth satellites**
  Earth satellite positions can be generated
  from public TLE elements describing their current orbit,
  which you can download from Celestrak. :doc:`Read more <earth-satellites>`

  .. TODO - update the following code with new approach

  .. testsetup::

    tle_text = """
    ISS (ZARYA)             
    1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082
    2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473
    """

  .. testcode::

    from skyfield.api import now, earth

    jd = now()
    boston = earth.topos('42.3583 N', '71.0603 W')
    satellite = earth.satellite(tle_text)

    # Geocentric

    apparent = satellite.gcrs(jd)

    # Topocentric

    apparent = boston.gcrs(jd).observe(satellite)

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

    # BCRS positions of Earth and Venus

    from skyfield.api import load

    planets = load('de421.bsp')

    print(planets('earth').at(utc=(1980, 1, 1)).position.au)
    print(planets('venus').at(utc=(1980, 1, 1)).position.au)

.. testoutput::

    [-0.16287311  0.88787399  0.38473904]
    [ 0.73245691 -0.01760889 -0.05403112]

The coordinates shown above are measured
using the Astronomical Unit (au),
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
the position of Mars in the night sky of Earth
by simply subtracting these two positions
to generate the vector difference between them.
But that would ignore the fact that light takes several minutes
to travel between Mars and the Earth.
The image of Mars in our sky
does not show us where it *is*, right now,
but where it *was* — several minutes ago —
when the light now reaching our eyes or instruments
actually left its surface.

Correcting for the light-travel time
does not simply fix a minor inconvenience,
but reflects a very deep physical reality.
Not only the light from Mars,
but *all* of its physical effects,
arrive no faster than the speed of light.
As Mars tugs at the Earth with its gravity,
we do not get pulled in the direction of the “real” Mars —
we get tugged in the direction of its time-delayed image
hanging in the sky above us!

So Skyfield offers a :meth:`~Position.observe()` method
that carefully backdates the position of another object
to determine where it was when it generated the image
that we see in our sky:

.. testcode::

    # Observing Mars from the Earth's position

    astrometric = planets('earth').observe('mars').at(utc=(1980, 1, 1))
    print(astrometric.position.au)

.. testoutput::

    [-0.92909581  0.21939949  0.15266885]

This light-delayed position is called the *astrometric* position,
and is traditionally mapped on a star chart
by the angles *right ascension* and *declination*
that you can compute using the :meth:`~Position.radec()` method
and display using their :meth:`~Angle.hstr()`
and :meth:`~Angle.dstr()` methods:

.. testcode::

    # Astrometric RA and declination

    ra, dec, distance = astrometric.radec()
    print(ra.hstr())
    print(dec.dstr())
    print(distance)

.. testoutput::

    11h 06m 51.22s
    +09deg 05' 09.2"
    0.96678 au

As we will explore in the next section,
objects never appear at exactly the position in the sky
predicted by the simple and ideal astrometric position.
But it is useful for mapping the planet
against the background of stars in a
`printed star atlas <http://www.amazon.com/s/?_encoding=UTF8&camp=1789&creative=390957&linkCode=ur2&pageMinusResults=1&suo=1389754954253&tag=letsdisthemat-20&url=search-alias%3Daps#/ref=nb_sb_noss_1?url=search-alias%3Daps&field-keywords=star%20atlas&sprefix=star+%2Caps&rh=i%3Aaps%2Ck%3Astar%20atlas&sepatfbtf=true&tc=1389754955568>`_,
because star atlases also use astrometric positions.

.. _apparent:

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
  for example, or is nearly behind Mars.
  The effect is small,
  but must be taken into account for research-grade results.

*Aberration*
  The velocity of the Earth itself through space
  adds a very slight slant to light arriving at our planet,
  in the same way that rain or snow
  seen through the windshield while driving
  appears to be slanting towards you
  because of your own motion.
  The effect is small enough — at most about 20 arcseconds —
  that it was not both observed and explained until 1728,
  when James Bradley realized that it provided the long-awaited proof
  that the Earth is indeed in motion in an orbit around the Sun.

Skyfield lets you apply both of these effects
by invoking the :meth:`~Astrometric.apparent()` method
on an astrometric position.
Like an astrometric position, an apparent position
is typically expressed as the angles
*right ascension* and *declination*:

.. testcode::

    # Apparent GCRS ("J2000.0") coordinates

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()

    print(ra.hstr())
    print(dec.dstr())
    print(distance)

.. testoutput::

    11h 06m 51.75s
    +09deg 05' 04.7"
    0.96678 au

But it is actually unusual to print apparent coordinates
in a permanent unchanging reference frame like the ICRS,
so you are unlikely to find the two values above
if you look up the position of Mars on 1980 January 1
in an almanac or by using other astronomy software.

Instead, apparent positions are usually expressed
relative to the Earth’s real equator and poles
as its rolls and tumbles through space —
which, after all,
is how right ascension and declination were defined
through most of human history,
before the invention of the ICRS axes.
The Earth’s equator and poles move at least slightly every day,
and move by very large amounts as years add up to centuries.

To ask for right ascension and declination
relative to the real equator and poles of Earth,
and not the ideal permanent axes of the ICRS,
simply add the keyword argument ``epoch='date'``
when you ask the apparent position for coordinates:

.. testcode::

    # Coordinates relative to true equator and equinox

    ra, dec, distance = apparent.radec(epoch='date')

    print(ra.hstr())
    print(dec.dstr())
    print(distance)

.. testoutput::

    11h 05m 48.68s
    +09deg 11' 35.7"
    0.96678 au

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
that we have been generating so far:

.. testcode::

    alt, az, distance = apparent.altaz()

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: to compute an apparent position, you must observe from a specific Earth location that you specify using a Topos instance

Instead, you have to give Skyfield your geographic location.
Astronomers use the term *topocentric*
for a position measured relative to a specific location on Earth,
so Skyfield represents Earth locations using :class:`Topos` objects
that you can generate by using the :meth:`Earth.topos` method
of an Earth object:

.. testcode::

    # Altitude and azimuth in the sky of a
    # specific geographic location

    boston = planets('earth').topos('42.3583 N', '71.0603 W')
    astro = boston.observe('mars').at(utc=(1980, 3, 1))
    app = astro.apparent()

    alt, az, distance = app.altaz()
    print(alt.dstr())
    print(az.dstr())
    print(distance)

.. testoutput::

    24deg 39' 49.8"
    93deg 13' 17.3"
    0.678874 au

So Mars was more than 24° above the horizon for Bostonians
on 1980 March 1 at midnight UTC.

The altitude returned from a plain :meth:`~Apparent.altaz()` call
is the ideal position
that you would observe if the Earth had no atmosphere.
You can also ask Skyfield to estimate
where an object might actually appear in the sky
after the Earth’s atmosphere has *refracted* their image higher.
If you know the weather conditions, you can specify them.

.. testcode::

    alt, az, distance = app.altaz(temperature_C=15.0,
                                  pressure_mbar=1005.0)
    print(alt.dstr())

.. testoutput::

    24deg 41' 55.8"

Or you can ask Skyfield to use a standard temperature and pressure
when generating its rough simulation of the effects of refraction.

.. testcode::

    alt, az, distance = app.altaz('standard')
    print(alt.dstr())

.. testoutput::

    24deg 41' 58.6"

Keep in mind that these are simply guesses.
The effects of the atmosphere,
with its many layers of heat and cold and wind and weather,
cannot be accurately modeled or predicted.
And note that refraction is only applied to objects above the horizon.
Objects below −1.0° altitude are not adjusted for refraction.
