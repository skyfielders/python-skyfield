
=============
 Coordinates
=============

Add text to Positions:

at end of radec
If instead you are starting with ra/dec coordinates
and want to convert them to a Skyfield position,
see _.

And do that with other ra/dec? and altaz?

“Coordinates” are the numbers we use to name astronomical positions.
For any given position,
Skyfield can return
any of at least a half-dozen different sets of coordinates
naming that position.

Through the discussion that follows,
remember that coordinates are merely names for a position.
Asking for a different coordinate doesn’t change the position itself.
As explained in the :doc:`positions` chapter,
there is a real physical difference between,
say, an astrometric position and an apparent position:
they are two distinct points in the sky.
But for a single position —
a single point in the sky —
all the various coordinates that Skyfield can produce
will identify the exact same spot,
and will point a telescope in its direction.

Cartesian coordinates versus Spherical coordinates
==================================================

A given coordinate system can be used to produce
two kinds of coordinate:
Cartesian and spherical.

*Cartesian* coordinates,
named after their inventor René Descartes,
measure a position as three distances |xyz|
along axes that stand at right angles to each other.

*Spherical* coordinates instead specify a location
using two angles plus a distance.
Just as geographic latitude and longitude
designate a position on the Earth,
the two angles of spherical coordinates
pick out a particular point on the celestial sphere —
a particular direction into the heavens —
and the distance says how far you would need to travel in that direction
to reach the target.

Through the following sections,
we will see how to ask for both Cartesian and spherical coordinates
in a variety of coordinate systems.

ICRS and equatorial right ascension and declination
===================================================

*Equatorial* coordinates,
which measure angles from the Earth’s equator and poles,
have already been described in the :doc:`positions` chapter:

* :ref:`positions:Right Ascension and Declination: Astrometric`
* :ref:`positions:Right Ascension and Declination: Apparent`

Here’s how the axes and angles are defined:

 | **Equatorial Coordinates**
 | *xy*-plane: Earth’s equator
 | *x*-axis: March equinox
 | *z*-axis: North celestial pole
 | ↔ Right ascension 0\ |h|–24\ |h| east around the equator
 | ↕ Declination ±90° from equator toward poles

There’s a problem
with the definition of this coordinate system:
thanks to `precession <https://en.wikipedia.org/wiki/Axial_precession>`_,
the Earth’s poles are in constant — if gradual — motion.
So equatorial coordinates
always have to specify
whether they are using the Earth’s equator and equinox
as of 1950, or 2000, or some other specific year.
In the late twentieth century
“J2000” coordinates became popular,
measured against where the Earth’s mean equator would be pointing
at noon Terrestrial Time on 2000 January 1.

But in 1998,
astronomers decided it was time to get off of the equatorial treadmill.
They defined a new coordinate system:
the ICRS (International Celestial Reference System).
Its axes point almost exactly along the old J2000 axes,
but are permanent and immobile and can achieve unbounded precision:
instead of moving with the Earth,
its axes are measured against the locations of very distance quasars.

In the :doc:`positions` chapter
you’ve already seen
how to retrieve ICRS coordinates in both Cartesian and spherical form:

.. testcode::

    # Let’s build an example position using
    # what we learned in the Positions chapter.

    from skyfield.api import load

    ts = load.timescale()
    t = ts.utc(2021, 12, 2, 14, 7)
    eph = load('de421.bsp')
    earth, mars = eph['earth'], eph['mars']
    position = earth.at(t).observe(mars)

    print('Cartesian ICRS:')

    x, y, z = position.xyz.au

    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    print('Spherical ICRS:')

    ra, dec, distance = position.radec()

    print(' ', ra, 'right ascension')
    print(' ', dec, 'declination')
    print(' ', distance, 'distance')

.. testoutput::

    Cartesian ICRS:
      x = -1.521 au
      y = -1.800 au
      z = -0.772 au

    Spherical ICRS:
      15h 19m 15.08s right ascension
      -18deg 08' 37.0" declination
      2.47983 au distance

These are the most popular coordinates
to use with astrometric positions.

But what if instead of using the ICRS
you want to pay attention to precession,
and get coordinates measured against the Earth’s true equator and equinox?
In that case,
retrieving spherical coordinates is easy:
simply provide a date to :meth:`~skyfield.positionlib.ICRF.radec()`
and it will return a right ascension and declination
relative to the Earth’s equator and equinox as of that date.

It’s less usual for someone to want |xyz| coordinates
measured against the real equator and equinox,
so there’s no built-in method to retrieve them.
Instead, we can import the Skyfield object
representing the reference frame
and pass it to the position’s
:meth:`~skyfield.positionlib.ICRF.frame_xyz()` method:

.. testcode::

    from skyfield.framelib import \
        true_equator_and_equinox_of_date as of_date

    print('Cartesian equinox-of-date coordinates:')

    x, y, z = position.frame_xyz(of_date).au

    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    print('Spherical equinox-of-date coordinates:')

    ra, dec, distance = position.radec(position.t)

    print(' ', ra, 'right ascension')
    print(' ', dec, 'declination')
    print(' ', distance, 'distance')

.. testoutput::

    Cartesian equinox-of-date coordinates:
      x = -1.510 au
      y = -1.808 au
      z = -0.775 au

    Spherical equinox-of-date coordinates:
      15h 20m 28.70s right ascension
      -18deg 13' 18.5" declination
      2.47983 au distance

Note that the distance is exactly the same as before,
because this is exactly the same position —
it’s merely being measured against a slightly different set of axes.

Horizonal coordinates
=====================

Altitude and azimuth have already been explained
in the :doc:`positions` chapter,
so you can start reading about them there:

* :ref:`positions:Azimuth and altitude from a geographic position`

The coordinate system is called *horizonal*
in the sense of “pertaining to the horizon.”

 | **Horizonal Coordinates**
 | *xy*-plane: Horizon
 | *x*-axis: North point on the horizon
 | *y*-axis: East point on the horizon (left-handed)
 | *z*-axis: Zenith
 | ↕ Altitude ±90° above or below horizon
 | ↔ Azimuth 0°–360° measured clockwise from north

As with the equatorial system,
the angles associated with horizontal coordinates are so popular
that Skyfield provides a built-in
method :meth:`~skyfield.positionlib.Apparent.altaz()` to retrieve them.
And you can generate |xyz| coordinates
by calling :meth:`~skyfield.positionlib.ICRF.frame_xyz()`
and passing it the geographic location itself as the reference frame:

.. testcode::

    # From the chapter on Positions:
    # computing altitude and azimuth.

    from skyfield.api import load, wgs84

    bluffton = wgs84.latlon(+40.8939, -83.8917)
    astrometric = (earth + bluffton).at(t).observe(mars)
    position = astrometric.apparent()

    print('Cartesian:')

    x, y, z = position.frame_xyz(bluffton).au

    print('  x = {:.3f} au north'.format(x))
    print('  y = {:.3f} au east'.format(y))
    print('  z = {:.3f} au up'.format(z))
    print()

    print('Spherical:')

    alt, az, distance = position.altaz()

    print('  Altitude:', alt)
    print('  Azimuth:', az)
    print('  Distance:', distance)

.. testoutput::

    Cartesian:
      x = -1.913 au north
      y = 1.200 au east
      z = 1.025 au up

    Spherical:
      Altitude: 24deg 24' 20.2"
      Azimuth: 147deg 54' 27.9"
      Distance: 2.47981 au

Note that some astronomers use the term “elevation”
for what Skyfield calls “altitude”:
the angle at which a target stands above the horizon.
Obviously both words are ambiguous,
since “elevation” can also mean a site’s vertical distance above sea level,
and since “altitude” can also mean an airplane’s height
above either sea level or the ground.

Hour Angle and Declination
==========================

If you are pointing a telescope or other instrument,
you might be interested in a variation on equatorial coordinates:
replacing right ascension with *hour angle,*
which measures ±180° from your own local meridian.

.. testcode::

    ha, dec, distance = position.hadec()

    print('Hour Angle:', ha)
    print('Declination:', dec, )
    print('Distance:', distance)

.. testoutput::

    Hour Angle: -02h 02m 28.94s
    Declination: -18deg 13' 16.4"
    Distance: 2.47981 au

To make the hour angle and declination even more useful
for pointing real-world instruments,
Skyfield includes the effect of polar motion
if you have :ref:`loaded a polar motion table <polar motion>`.
In that case the declination you get from
:meth:`~skyfield.positionlib.ICRF.hadec()`
will vary slightly from the declination returned by
:meth:`~skyfield.positionlib.ICRF.radec()`,
which doesn’t include polar motion.

ECI versus ECEF coordinates
===========================

Here’s a quick explanation of two acronyms
that you’re likely to run across in discussions about coordinates.

ECI stands for *Earth-Centered Inertial*
and specifies coordinates that are
(a) measured from the Earth’s center
and (b) that don’t rotate with the Earth itself.
The very first coordinates we computed in this chapter,
for example,
qualify as ECI coordinates,
because the ``position`` used the Earth as its center
and because the ICRS system of right ascension and declination
stays fixed on the celestial sphere
even as the Earth rotates beneath it.

ECEF stands for *Earth-Centered Earth-Fixed*
and specifies coordinates that are
(a) measured from the Earth’s center
but (b) which rotate with the Earth instead of staying fixed in space.
A fixed latitude and longitude on the Earth’s surface is a good example.
We will learn about generating ECEF coordinates in the next section.

Geographic ITRS latitude and longitude
======================================

Skyfield uses the standard ITRS reference frame
for specifying Earth-fixed positions
that are measured from the rotating Earth’s surface.

 | **ITRS Coordinates**
 | *xy*-plane: Earth’s equator
 | *x*-axis: 0° longitude on the equator
 | *y*-axis: 90° east longitude on the equator
 | *z*-axis: North pole
 | ↕ Latitude ±90° from equator toward poles
 | ↔ Longitude ±180° from prime meridian with east positive

The definition of latitude
depends on whether you model the Earth as a simple sphere
or more realistically as a slightly flattened ellipsoid.
The most popular choice today is to use the WGS84 ellipsoid,
which is the one used by the GPS system.

.. testcode::

    from skyfield.api import wgs84
    from skyfield.framelib import itrs

    # Important: must start with a position
    # measured from the Earth’s center.
    position = earth.at(t).observe(mars)

    print('Cartesian:')

    x, y, z = position.frame_xyz(itrs).au

    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    print('Geographic:')

    lat, lon = wgs84.latlon_of(position)
    height = wgs84.height_of(position)

    print(' {:.4f}° latitude'.format(lat.degrees))
    print(' {:.4f}° longitude'.format(lon.degrees))
    print(' {:.0f} km above sea level'.format(distance.km))

.. testoutput::

    Cartesian:
      x = 1.409 au
      y = -1.888 au
      z = -0.775 au

    Geographic:
     -18.2218° latitude
     -53.2660° longitude
     370974969 km above sea level

Note that height is measured from sea level,
not distance from the center of the Earth.

The code above is slightly inefficient,
because :meth:`~skyfield.toposlib.Geoid.height_of()`
will wind up recomputing several values
that were already computed in :meth:`~skyfield.toposlib.Geoid.latlon_of()`.
If you need both, it’s more efficient to call
:meth:`~skyfield.toposlib.Geoid.geographic_position_of()`.

There’s also a :meth:`~skyfield.toposlib.Geoid.subpoint_of()` method
if you want Skyfield to compute the geographic position
of the sea-level point beneath a given celestial object.

.. Once fully supported, illustrate round-trips like

    xyz = m.frame_xyz(itrs)
    from skyfield.positionlib import ICRS
    position = ICRS.from_time_and_frame_vectors(t, itrs, xyz, None)

Ecliptic coordinates
====================

*Ecliptic coordinates* are measured from the plane of the Earth’s orbit.
They are useful
when making maps and diagrams of the Solar System
and when exploring the properties of orbits around the Sun,
because it places the orbits of the major planets
nearly flat against the *xy*-plane —
unlike right ascension and declination,
which twists the Solar System up at the 23° tilt of the Earth’s own axis.

You might be tempted to ask
why we measure against the plane of the Earth’s orbit,
instead of averaging together all the planets
to compute the “invariable plane” of the whole Solar System
(to which the Earth’s orbit is inclined by something like 1.57°).
The answer is: precision.
We know the plane of the Earth’s orbit to many decimal places,
because the Earth carries all of our highest-precision observatories
along with it as it revolves around the Sun.
Our estimate of the invariable plane, by contrast,
is a mere average
that changes — at least slightly —
every time we discover a new trans-Neptunian object, asteroid, or comet.
So the Earth’s own orbit remains the best basis
for a coordinate system oriented to the Solar System.

 | **Ecliptic Coordinates**
 | *xy*-plane: Ecliptic plane
 | *x*-axis: March equinox
 | *z*-axis: North ecliptic pole
 | ↕ Latitude ±90° above or below the ecliptic
 | ↔ Longitude 0°–360° measured east from March equinox

.. testcode::

    from skyfield.framelib import ecliptic_frame

    print('Cartesian ecliptic coordinates:')

    x, y, z = position.frame_xyz(ecliptic_frame).au

    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    print('Spherical ecliptic coordinates:')

    lat, lon, distance = position.frame_latlon(ecliptic_frame)

    print(' {:.4f} latitude'.format(lat.degrees))
    print(' {:.4f} longitude'.format(lon.degrees))
    print(' {:.3f} au distant'.format(distance.au))

.. testoutput::

    Cartesian ecliptic coordinates:
      x = -1.510 au
      y = -1.967 au
      z = 0.007 au

    Spherical ecliptic coordinates:
     0.1732 latitude
     232.4801 longitude
     2.480 au distant

Note the very small values returned
for the ecliptic *z* coordinate
and for the ecliptic latitude —
because we are measuring against the plane
in which both the Earth and Mars revolve around the Sun.

Galactic coordinates
====================

*Galactic coordinates* are measured
against the plane and center of our Milky Way galaxy —
or at least as best as we can approximate the galaxy’s sprawling structure
from our vantage point here deep inside the Orion Arm:

 | **Galactic Coordinates**
 | *xy*-plane: Galactic plane
 | *x*-axis: Galactic center
 | *z*-axis: North galactic pole
 | ↕ Latitude ±90° above galactic plane
 | ↔ Longitude 0°–360° east from galactic center

.. testcode::

    from skyfield.framelib import galactic_frame

    print('Cartesian galactic coordinates:')

    x, y, z = position.frame_xyz(galactic_frame).au

    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    print('Spherical galactic coordinates:')

    lat, lon, distance = position.frame_latlon(galactic_frame)

    print(' {:.4f} latitude'.format(lat.degrees))
    print(' {:.4f} longitude'.format(lon.degrees))
    print(' {:.3f} au distant'.format(distance.au))

.. testoutput::

    Cartesian galactic coordinates:
      x = 2.029 au
      y = -0.527 au
      z = 1.324 au

    Spherical galactic coordinates:
     32.2664 latitude
     345.4330 longitude
     2.480 au distant

Astronomers have generated a series of more and more precise estimates
of our galaxy’s orientation over the past hundred years.
Skyfield uses the `IAU 1958 Galactic System II
<https://adsabs.harvard.edu/full/1960MNRAS.121..123B>`_,
which is believed to be accurate to within ±0.1°.

.. TODO section on Velocity

Turning coordinates into a position
===================================

.. TODO

All of the above examples take a Skyfield position and return coordinates,
but sometimes you start with coordinates in a reference frame
and want to produce a position.
If you happen to start with |xyz| coordinates,
then you can create a position with a method call:

    #ICRS
    from skyfield.framelib import ecliptic
    from skyfield.positionlib import Apparent

    t = ?
    xyz = [...]
    a = Apparent.from_time_and_frame_vectors(t, ecliptic, xyz, None)

If instead of a Cartesian |xyz| vector
you start with a right ascension and declination,
then 

If your coordinates are expressed as some other pair of angles,
then you will have to dig a bit deeper into Skyfield
/and do the conversion yourself/
for the routine to convert between the two.

    from skyfield.functions import from_spherical

    lat, lon, distance =
    xyz = from_spherical(distance.au, lat.radians, lon.radians)

Then, use the same maneuver shown above
to turn the |xyz| vector into a Skyfield position.


position_of_radec
from_time_and_frame_vectors
