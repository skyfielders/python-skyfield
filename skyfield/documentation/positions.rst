
===========
 Positions
===========

.. currentmodule:: skyfield.positionlib

Skyfield stores the locations of celestial bodies
as |xyz| vectors that it calls *positions.*
There are several kinds of celestial object
for which Skyfield can produce a position.
If there’s a particular kind of object you’re interested in,
you might want to read its documentation first,
then return here when you’re ready
to learn more about the position objects that they generate:

* `planets`
* `stars`
* `earth-satellites`
* `kepler-orbits` (comets and asteroids)

You can also build a position object yourself
by providing |xyz| coordinates to a position class:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup()

.. TODO add time to this example?

.. testcode::

    from skyfield.positionlib import Barycentric

    x = 3141.0
    y = 2718.0
    z = 5820.0
    barycentric = Barycentric([x, y, z])

This document focuses
on what you can do with positions once they’ve been computed,
and on how position objects are used in Skyfield code.

Quick reference
===============

For a succinct catalog of the attributes and methods
offered by Skyfield’s position classes,
follow these links to the API documentation.

* :ref:`api:Astronomical positions`

  All the attributes and methods of the core Skyfield position classes.

* :ref:`api:Units`

  Skyfield’s distance, velocity, and angle classes,
  which offer simple attributes like ``km``
  but also fancy methods like ``dstr()``
  that formats an angle as degrees, minutes, and seconds.

.. _ICRS:

The ICRS reference system and J2000
===================================

.. TODO Make ra/dec and alt/az links to subsequent sections

Even though Skyfield scripts often produce output
in spherical coordinates —
like right ascension and declination,
or altitude and azimuth —
Skyfield always stores positions internally
as Cartesian |xyz| vectors
oriented along the axes of the International Celestial Reference System (ICRS).

The ICRS is a higher-accuracy replacement
for the old J2000 reference system.
It’s defined using the positions of very distant quasars,
so its precision can improve each decade
as radio telescopes measure quasar positions ever more accurately.

If you want to interact directly with ICRS |xyz| coordinates,
here is where its axes point:

* *x-axis* — Aims at 2000 January 1 position of the Vernal Equinox,
  which is defined more technically
  as the ascending node of the ecliptic on the mean celestial equator.
  Ancient astronomers called this “the first point of Ares”
  but precession has gradually shifted it into the constellation Pisces.

* *y-axis* — Aims at the point 90° east of the Vernal Equinox
  along the celestial equator,
  which lies south of Betelgeuse and a few degrees east of Orion’s belt.

* *z-axis* — Aims at the North Celestial Pole.

The ICRS axes are within 0.02 arcseconds of the old J2000 axes,
so many scripts simply treat J2000 coordinates
as modern ICRS coordinates.

Barycentric → Astrometric → Apparent
====================================

The most common calculations in Skyfield
produce a succession of three |xyz| positions,
one for the observer and two for the body being observed.
The three positions can be hard to find in Skyfield code
because they are often created and used
without ever being given a name.
As an example,
let’s use an :doc:`ephemeris <planets>`
to compute a geocentric position for Mars:

.. testcode::

    from skyfield.api import load

    ts = load.timescale()
    t = ts.utc(1980, 1, 1)
    planets = load('de421.bsp')
    earth, mars = planets['earth'], planets['mars']

    # Three positions in a single line of code!
    d = earth.at(t).observe(mars).apparent().distance()

    print('Mars is {:.2f} au from Earth'.format(d.au))

.. testoutput::

    Mars is 0.97 au from Earth

By chaining together four different methods,
the line that computes the distance ``d = …``
creates and discards three different Skyfield positions
in a single line of code!
To learn about them,
the beginner will probably need to slow up and re-write that line
so each object is assigned a separate name:

.. testcode::

    # The “d =” line from the previous example,
    # rewritten to give each position a name.

    barycentric = earth.at(t)
    astrometric = barycentric.observe(mars)
    apparent = astrometric.apparent()
    d = apparent.distance()

This is a common Python pattern.
By assigning names to intermediate values,
the programmer
— without changing the code’s result —
can pivot between
succinct code that fits on a single line
and more verbose code that names each intermediate value.

Now that we’ve given them names,
we can discuss the three positions:

* A :class:`Barycentric` position
  measures from the Solar System’s center of mass.
  This places its |xyz| vector
  in the Barycentric Celestial Reference System (BCRS) —
  a frame of reference that’s inertial enough
  to support the :meth:`~Barycentric.observe()` method.

  You’ll usually start a Skyfield script
  by generating a barycentric position
  for the *center* from which you’ll be observing —
  whether that’s the Earth,
  or a specific location on the Earth’s surface,
  or another body like a satellite, planet, or moon.

* An :class:`Astrometric` position
  is returned by the :meth:`Barycentric.observe()` method which,
  given a target you want to observe,
  applies the effect of light travel time.
  For example,
  on Earth we see the Moon where it was about 1.3 seconds ago,
  the Sun where it was 8 minutes ago,
  Jupiter where it was more than half an hour ago,
  and Neptune where it was about 4 hours ago.

.. TODO put this below
   The word “astrometric” is Greek for “star measure”
   because it’s the position
   where you would draw the target body on a star chart.

* An :class:`Apparent` position is computed
  by calling the :meth:`Astrometric.apparent()` method.
  This applies two real-world effects
  that slightly shift everything in the sky:
  the aberration of light
  produced by the observer’s own motion through space,
  and the gravitational deflection of light
  that passes close to masses like the Sun and Jupiter —
  and, for an observer on the Earth’s surface,
  for deflection produced by the Earth’s own gravity.
  The result is an “apparent” position
  telling you where the target will really “appear” in tonight’s sky;
  it’s the direction you should point your telescope.

  When an apparent position is measured from the Earth’s center,
  it can be described more formally
  as an |xyz| position
  in Geocentric Celestial Reference System (GCRS) coordinates.

Position attributes
===================

Five basic attributes are available on each position:

.. include:: api.rst
   :start-after: PAT START
   :end-before: PAT END

The first three attributes listed above
are simple instances of Skyfield’s distance, velocity, and time classes,
which you can learn more about
by clicking on their class names.
They support operations like:

.. testcode::

    print('Earth x,y,z:', barycentric.position.au, 'au')
    print('Mars relative velocity:', astrometric.velocity.km_per_s, 'km/s')
    print('Time of observation:', apparent.t.utc_strftime())

.. testoutput::

    Earth x,y,z: [-0.16287311  0.88787399  0.38473904] au
    Mars relative velocity: [12.66638873 -8.12551301 -3.38356109] km/s
    Time of observation: 1980-01-01 00:00:00 UTC

Note that the distance unit attributes
like ``au`` and ``km``
and the velocity unit attributes
like ``au_per_d`` and ``km_per_s``
are each three-element NumPy arrays offering |xyz| coordinates.
If the times you’re using are :ref:`date-arrays`,
then the distance and velocity will each have an additional dimension
offering as many |xyz| coordinates as there are dates in your array.

The ``.center`` and ``.target`` attributes require a bit more explanation.
They specify the origin and destination of the vector,
and are displayed if you ask Python to print the vector:

.. testcode::

    print(apparent)

.. testoutput::

    <Apparent GCRS position and velocity at date t center=399 target=499>

.. TODO document NAIF codes and have brief summary, on ephemeris page

In this case our positions were generated from a JPL ephemeris file,
so the center and target are the simple integers
``399`` indicating the Earth and ``499`` indicating Mars.
If instead the center or target were defined by a Skyfield object
like an Earth satellite or latitude-longitude position,
then the ``.center`` and ``.target``
will be those Skyfield objects themselves.

Right ascension and declination: astrometric
============================================

The most popular coordinate system for star catalogs
treats the night sky as a slowly spinning globe
that we view from the inside.
Just as we specify position on the Earth’s globe
by degrees longitude that you would travel along the equator
and then degrees latitude that you would travel north or south,
coordinates in the sky are measured by two angles.
These are called “equatorial coordinates”
since they are measured with respect to the equator.

* *Right ascension* (“RA”)
  is the sky’s equivalent of longitude,
  and is measured east along the celestial equator.
  Since the sky makes roughly one complete turn every day,
  RA is usually expressed in units of 24 hours
  rather than 360 degrees.
  This supports casual inferences,
  like the fact
  that a star with an RA of 3\ |h| will climb to your meridian
  one hour later than a star with an RA of 2\ |h|.

  This coordinate presented astronomers
  with the same challenge that longitude presented geographers:
  the arbitrary choice of a starting point.
  In the case of longitude,
  we measure east and west
  from the Airy Transit Circle at the Royal Observatory Greenwich.
  For RA,
  astronomers measure east from the Vernal Equinox,
  the point where the Sun crosses the celestial equator in March
  as it passes from the southern to the northern half of the sky.

* *Declination* (“Dec”)
  is the sky’s equivalent of latitude,
  measured north and south of the celestial equator
  in degrees, with north being positive.
  The North Celestial Pole is at +90°
  and the South Celestial Pole is at −90°.

There are two common uses for RA/Dec coordinates.

The first common use is as a reference coordinate system
for perennial knowledge like star charts,
celestial catalogs, and record keeping.
This use case presents two complications.

The first complication
is that real star positions appear to be in constant motion.
The Earth’s revolution around the Sun
moves the stars in little circles thanks to the aberration of light,
and their light can be bent and their images shifted
by the gravity of the Sun, Jupiter, and even the Earth.
Because these effects are periodic or temporary,
they are entirely unsuitable for star charts,
which are supposed to be independent of any particular day of the year.

For this reason,
always generate reference RA/Dec coordinates
from astrometric Skyfield positions,
never from apparent positions.

The second complication is that the Earth’s poles
don’t point in a fixed direction
but gradually trace 26,000-year circles around the sky.
This requires RA/Dec coordinates to specify which year’s “equinox”
their right ascension is measured from,
which will also be the year whose poles and equator they use.

The modern standard for astrometric RA/Dec is the ICRS
(:ref:`described above <ICRS>`).
With axes that are fixed in the directions of the J2000 equinox and poles,
the ICRS is a permanent coordinate system
that will never suffer precession.
Skyfield returns ICRS coordinates
if you simply call :meth:`~ICRF.radec()`
without an argument:

.. testcode::

    # Astrometric RA/Dec.
    ra, dec, distance = astrometric.radec()

    print('RA:', ra)
    print('Dec:', dec)
    print('Distance:', distance)

.. testoutput::

    RA: 11h 06m 51.22s
    Dec: +09deg 05' 09.2"
    Distance: 0.96678 au

.. TODO add docstring to Angle as gentle introduction for API docs

Here we have printed RA and declination using their default units.
See the :class:`~skyfield.units.Angle` class documentation
for the attributes and methods
by which you can retrieve and format an angle’s value
in units of your own choice.

If your project specifically requires coordinates
expressed in the RA/Dec of an older equinox,
you can build a time object and pass it to  :meth:`~ICRF.radec()`:

.. testcode::

    # Astrometric RA/Dec relative to another equinox:
    # the J1991.25 epoch used by the Hipparcos catalog.

    equinox = ts.J(1991.25)
    ra, dec, distance = astrometric.radec(equinox)

Right ascension and declination: apparent
=========================================

The other reason that you might generate RA/Dec is practical:
you are planning to point a telescope,
and after orienting its equatorial mount to the Earth’s poles
you want a declination that stays fixed all night
while your telescope tracks the specified right ascension across the sky.

In this case you aren’t likely to be satisfied
with RA/Dec coordinates of some other era.
You’ll want them measured against where the Earth’s poles
are really pointing tonight.

Skyfield uses the high precision IAU 2000A model
to compute both the precession that carries the Earth’s poles
in their 26,000-year circle around the sky
and also the short term wobbles in the Earth’s orientation
that are called nutation.
Together they give Skyfield an accurate assessment
of the direction of the poles and equator
on a given date.

When pointing a telescope,
always use apparent coordinates,
since you will want every possible real-world effect accounted for.

While you could pass each position’s time
to its own :meth:`~ICRF.radec()` method,
that would be a bit tedious,
so Skyfield provides a shortcut:
if you pass the string ``'date'``
then the RA/Dec coordinates use the equinox and poles
of the date of the position itself:

.. testcode::

    # Apparent RA/Dec.
    ra, dec, distance = apparent.radec('date')

    print('RA:', ra)
    print('Dec:', dec)
    print('Distance:', distance)

.. testoutput::

    RA: 11h 05m 48.68s
    Dec: +09deg 11' 35.7"
    Distance: 0.96678 au

.. “is compounded of?”

Note that there are two sources of difference
between the astrometric coordinates printed in the previous section
and the apparent coordinates printed here.
We switched to the RA/Dec system of a different year,
so even the exact same position
would have been assigned a different coordinate than before.
But we have also asked for the RA/Dec of a whole different point in the sky:
the point where Mars will actually appear,
not the ideal point where it should sit on a star chart
that ignores aberration and deflection.

.. TODO can generate observer position only once and use for N .observe() calls

.. TODO move this to an accuracy page entry on aberration

  The velocity of the Earth itself through space
  adds a very slight slant to light arriving at our planet,
  in the same way that rain or snow
  seen through the windshield while driving
  appears to be slanting towards you
  because of your own motion.
  The effect is small enough — at most about 20 arcseconds —
  that only in 1728 was it finally observed and explained,
  when James Bradley realized that it provided the long-awaited proof
  that the Earth is indeed in motion in an orbit around the Sun.

Azimuth and altitude from a geographic position
===============================================

The final result that many users seek
is the altitude and azimuth of an object
above their own local horizon.

* *Altitude* measures the angle above or below the horizon.
  The zenith is at +90°,
  an object on the horizon’s great circle is at 0°,
  and the nadir beneath your feet is at −90°.

* *Azimuth* measures the angle around the sky from the north pole:
  0° means exactly north, 90° is east, 180° is south, and 270° is west.

Altitude and azimuth are computed
by calling the :meth:`~Apparent.altaz()` method on an apparent position.
But because the method needs to know which local horizon to use,
it does not work
on the plain geocentric (“Earth centered”) positions
that we have been generating so far:

.. testcode::

    alt, az, distance = apparent.altaz()

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: to compute an apparent position, you must observe from a specific Earth location that you specify using a latitude and longitude

You can specify a location on Earth
by giving its latitude and longitude
to a standard “geodetic system” that models the Earth’s shape.
The most popular model is WGS84,
used by most modern maps
and also by the Global Positioning System (GPS).
If you are given a longitude and latitude without a datum specified,
they are probably WGS84 coordinates.

You can pass the latitude and longitude to a datum like WGS84
by calling its :meth:`~skyfield.toposlib.Geoid.latlon()` method.
Skyfield provides constants ``N``, ``S``, ``E``, and ``W``
that are each positive one or negative one,
in case you don’t want to remember which directions are positive.

.. testcode::

    from skyfield.api import N,S,E,W, wgs84

    # Altitude and azimuth in the sky of a
    # specific geographic location

    boston = earth + wgs84.latlon(42.3583 * N, 71.0603 * W, elevation_m=43)
    astro = boston.at(ts.utc(1980, 3, 1)).observe(mars)
    app = astro.apparent()

    alt, az, distance = app.altaz()
    print(alt.dstr())
    print(az.dstr())
    print(distance)

.. testoutput::

    24deg 30' 27.2"
    93deg 04' 29.5"
    0.678874 au

So Mars was more than 24° above the horizon for Bostonians
on 1980 March 1 at midnight UTC.

The altitude returned from a plain :meth:`~Apparent.altaz()` call
is the ideal position
that you would observe if the Earth had no atmosphere.
You can also ask Skyfield to estimate
where an object might actually appear in the sky
after the Earth’s atmosphere has *refracted* its image higher.
If you know the weather conditions, you can specify them.

.. testcode::

    alt, az, distance = app.altaz(temperature_C=15.0,
                                  pressure_mbar=1005.0)
    print(alt.dstr())

.. testoutput::

    24deg 32' 34.1"

Or you can ask Skyfield to use a standard temperature and pressure
when generating its rough simulation of the effects of refraction.

.. testcode::

    alt, az, distance = app.altaz('standard')
    print(alt.dstr())

.. testoutput::

    24deg 32' 36.4"

Keep in mind
that the computed effect of refraction is simply an estimate.
The effects of your local atmosphere,
with its many layers of heat and cold and wind and weather,
cannot be predicted to high precision.
And note that refraction is only applied to objects above the horizon.
Objects below −1.0° altitude are not adjusted for refraction.

Comparing positions
===================

If you want to know the angle between two positions in the sky,
call the
:meth:`~skyfield.positionlib.ICRF.separation_from()`
method of one of the positions
and pass it the other position as its argument.
The result will be an :class:`~skyfield.units.Angle` object.

If instead you want to know the distance between two positions,
subtract the position you want to use as the starting point
from the other position.
The result of this vector math will also itself be a position vector:

.. testcode::

    v = planets['moon'].at(t) - planets['earth'].at(t)
    print('The Moon is %d km away' % v.distance().km)

.. testoutput::

    The Moon is 385010 km away

Two vectors do not need to have the same date and time to be subtracted.
By comparing two positions with different times,
you can measure how far an object has moved:

.. testcode::

    t1 = ts.utc(2015, 10, 11, 10, 30)
    t2 = ts.utc(2015, 10, 11, 10, 31)

    moon = planets['moon']
    p1 = moon.at(t1)
    p2 = moon.at(t2)

    km = (p2 - p1).distance().km
    print('In one minute the Moon moved %d km' % km)

.. testoutput::

    In one minute the Moon moved 1736 km

.. _reference_frames:

Coordinates in other reference frames
=====================================

You can ask Skyfield to express a position
in a number of other reference frames
besides the standard ICRF reference frame
(the modern equivalent to J2000 coordinates)
that Skyfield uses internally.
For example,
to express the position of the Moon relative to the rotating Earth:

.. testcode::

    from skyfield.framelib import itrs

    a = earth.at(t).observe(planets['moon']).apparent()
    x = a.frame_xyz(itrs)
    print(x.km.astype(int))

.. testoutput::

    [ 349045 -106774  122510]

The three position methods that accept a reference frame argument are:

* `frame_xyz()`
* `frame_xyz_and_velocity()`
* `frame_latlon()`

Position classes also support a constructor
that accepts a position vector and a velocity vector
in a particular reference frame:

* `from_time_and_frame_vectors()`

Here are the reference frames defined in the ``framelib`` module
(click on their names for more detailed descriptions):

* `true_equator_and_equinox_of_date`
* `itrs`
* :data:`~skyfield.framelib.ecliptic_J2000_frame`
* :data:`~skyfield.framelib.ecliptic_frame`
* :data:`~skyfield.framelib.galactic_frame`

See also :doc:`planetary` for reference frames
that are not included with Skyfield
but that you can load from NASA reference files.

.. testcleanup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()
