
=============
 Coordinates
=============

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
But given a single position —
a single point in the sky —
all the various coordinates you can produce
mean exactly the same thing:
the direction you would point a telescope to put that point at its center.

Cartesian coordinates versus Spherical coordinates
==================================================

A given coordinate system can be used to produce
two kinds of coordinate:
Cartesian and spherical.

*Cartesian* coordinates
specify a position by providing three distances |xyz|
that are measured along axes that stand at right angles to each other.

*Spherical* coordinates instead specify a location
using two angles plus a distance.
Just as a geographic latitude and longitude
designate a position on the Earth,
so the two angles of a spherical coordinate
pick out a particular point in the sky —
a particular point on the celestial sphere —
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

* :ref:`positions:Astrometric right ascension and declination`
* :ref:`positions:Apparent right ascension and declination`

Here’s how the axes and angles are defined:

 | **Equatorial Coordinates**
 | *xy*-plane: Earth’s equator
 | *x*-axis: March equinox
 | *z*-axis: North celestial pole
 | ↔ Right ascension 0\ |h|–24\ |h| east around the equator
 | ↕ Declination ±90° from equator toward poles

There’s a problem with this coordinate system:
it slowly changes.
Thanks to `precession <https://en.wikipedia.org/wiki/Axial_precession>`_,
the Earth’s poles are in constant — if gradual — motion.
So equatorial coordinates
always have to specify
whether they are measuring against the Earth’s equator and equinox
as of 1950, or 2000, or some other specific year.
In the late twentieth century
“J2000” coordinates became popular,
measured against the Earth’s mean equator
as of noon Terrestrial Time on 2000 January 1.

But in 1998,
astronomers decided it was time to get off of the equatorial treadmill.
They defined a new coordinate system:
the ICRS (International Celestial Reference System).
Its axes point almost exactly along the old J2000 axes,
but are permanent and immobile,
and can achieve unbounded precision
because instead of using the Earth
its axes are tied to the locations of very distant quasars.

In the :doc:`positions` chapter
you’ve already seen
how to retrieve ICRS coordinates in both Cartesian and spherical form:

.. testcode::

    # Building an example position.

    from skyfield.api import load

    ts = load.timescale()
    t = ts.utc(2021, 12, 2, 14, 7)
    eph = load('de421.bsp')
    earth, mars = eph['earth'], eph['mars']
    position = earth.at(t).observe(mars)

    # Producing coordinates.

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

Note that ``.xyz`` is a :class:`~skyfield.units.Distance` object
and supports not only ``.au``
but also several other units of measure.

If instead of using the permanent ICRS
you want to measure coordinates against the Earth’s axes as they precess,
then retrieving spherical coordinates is easy:
simply provide a date to :meth:`~skyfield.positionlib.ICRF.radec()`
and it will return a right ascension and declination
relative to the Earth’s equator and equinox as of that date.

It’s less usual for someone to want |xyz| coordinates
measured against the real equator and equinox,
so there’s no built-in method to retrieve them.
Instead, you’ll need to import the reference frame itself
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

.. _horizontal-coordinates:

Altitude and azimuth (‘horizonal’ coordinates)
==============================================

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
method :meth:`~skyfield.positionlib.ICRF.altaz()` to retrieve them,
while |xyz| coordinates require a call to
:meth:`~skyfield.positionlib.ICRF.frame_xyz()`
with the geographic location itself passed as the reference frame:

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
      Altitude: 24deg 24' 20.6"
      Azimuth: 147deg 54' 28.8"
      Distance: 2.47981 au

Note that some astronomers use the term “elevation”
for what Skyfield calls “altitude”:
the angle at which a target stands above the horizon.
Obviously both words are ambiguous,
since “elevation” can also mean a site’s vertical distance above sea level,
and “altitude” can mean an airplane’s height
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

    Hour Angle: -02h 02m 28.88s
    Declination: -18deg 13' 16.4"
    Distance: 2.47981 au

To make the hour angle and declination even more useful
for pointing real-world instruments,
Skyfield includes the effect of polar motion
if you have :ref:`loaded a polar motion table <polar-motion>`.
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
An example would be the latitude and longitude of the Lowell Observatory,
which stays in place on the Earth’s surface
but from the point of view of the rest of the Solar System
is rotating with the Earth.
A fixed  on the Earth’s surface is a good example.
We will learn about generating ECEF coordinates in the next section.

Geographic ITRS latitude and longitude
======================================

Skyfield uses the standard ITRS reference frame
to specify positions
that are fixed relative to the Earth’s surface.

 | **ITRS Coordinates**
 | *xy*-plane: Earth’s equator
 | *x*-axis: 0° longitude on the equator
 | *y*-axis: 90° east longitude on the equator
 | *z*-axis: North pole
 | ↕ Latitude ±90° from equator toward poles
 | ↔ Longitude ±180° from prime meridian with east positive

A location’s latitude will vary slightly
depending on whether you model the Earth as a simple sphere
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
     -53.2662° longitude
     370974969 km above sea level

Note that height is measured from sea level,
not from the center of the Earth.

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

.. _ecliptic-coordinates:

Ecliptic coordinates
====================

*Ecliptic coordinates* are measured from the plane of the Earth’s orbit.
They are useful
when making maps and diagrams of the Solar System
and when exploring the properties of orbits around the Sun,
because they place the orbits of the major planets
nearly flat against the *xy*-plane —
unlike right ascension and declination,
which twist the Solar System up at a 23° angle
because of the tilt of the Earth’s axis.

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
every time we discover a new asteroid, comet, or trans-Neptunian object.
So the Earth’s own orbit winds up being the most pragmatic choice
for a coordinate system oriented to the Solar System.

 | **Ecliptic Coordinates**
 | *xy*-plane: Ecliptic plane (plane of Earth’s orbit)
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
and for the ecliptic latitude,
because Mars revolves around the Sun
in very nearly the same plane as the Earth.

Galactic coordinates
====================

*Galactic coordinates* are measured
against the disc of our own Milky Way galaxy,
as measured from our vantage point here inside the Orion Arm:

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

Velocity
========

If you need not only a position vector but a velocity vector
relative to a particular reference frame,
then switch from
the :meth:`~skyfield.positionlib.ICRF.frame_xyz()` method
to the :meth:`~skyfield.positionlib.ICRF.frame_xyz_and_velocity()` method.
As explained in its documentation,
its return values include the velocity vector.

Turning coordinates into a position
===================================

All of the above examples take a Skyfield position and return coordinates,
but sometimes you start with coordinates
and want to produce a position.

If you happen to start with ICRS |xyz| coordinates,
then you can create a position with a function call:

.. testcode::

    from skyfield.positionlib import build_position

    icrs_xyz_au = [-1.521, -1.800, -0.772]
    position = build_position(icrs_xyz_au, t=t)

But it’s probably more common for you
to have been given a target’s right ascension and declination.
One common approach is to create a
:class:`~skyfield.starlib.Star`
as described in the :doc:`stars` chapter,
which will give you an object
that you can pass to ``.observe()``
like any other Skyfield body.
But you can also create a position directly
by using the ``from_radec()`` method carried by each position class.
To create an apparent position, for example:

.. testcode::

    from skyfield.positionlib import Apparent

    position = Apparent.from_radec(ra_hours=5.59, dec_degrees=5.45)

A final common situation
is that you measured an altitude and azimuth relative to your horizon
and want to learn the right ascension and declination of that position.
The solution is to use the :meth:`~skyfield.positionlib.ICRF.from_altaz()`
method,
but there’s a catch:
because the true coordinates behind any particular altitude and azimuth
are changing every moment as the Earth spins,
you first need to compute the position of your observatory
``.at()`` the moment you measured the altitude and azimuth,
and only then call the method.

.. testcode::

    # What are the coordinates of the zenith?

    b = bluffton.at(t)
    apparent = b.from_altaz(alt_degrees=90.0, az_degrees=0.0)

    ra, dec, distance = apparent.radec()
    print('Right ascension:', ra)
    print('Declination:', dec)

.. testoutput::

    Right ascension: 13h 17m 00.26s
    Declination: +41deg 00' 27.7"

If you find yourself in an even less common situation,
like needing to build a position from ecliptic or galactic coordinates,
then —
while there aren’t yet any documented examples for you to follow —
you might be able to assemble a solution together from these pieces:

* The position constructor method
  :meth:`~skyfield.positionlib.ICRF.from_time_and_frame_vectors()`.
* The ``from_spherical(r, theta, phi)`` method in ``skyfield/functions.py``.

.. TODO

 If your coordinates are expressed as some other pair of angles,
 then you will have to dig a bit deeper into Skyfield
 /and do the conversion yourself/
 for the routine to convert between the two.

    from skyfield.functions import from_spherical

    lat, lon, distance =
    xyz = from_spherical(distance.au, lat.radians, lon.radians)

 Then, use the same maneuver shown above
 to turn the |xyz| vector into a Skyfield position.

 from_time_and_frame_vectors

Rotation Matrices
=================

If you are doing some of your own mathematics,
you might want access to the low-level 3×3 rotation matrices
that define the relationship between each coordinate reference frame
and the ICRS.
To compute a rotation matrix,
simply pass a time to the frame’s ``rotation_at()`` method:

.. testcode::

    # 3×3 rotation matrix: ICRS → frame

    R = ecliptic_frame.rotation_at(t)
    print(R)

.. testoutput::

    [[ 0.99998613 -0.00482998 -0.00209857]
     [ 0.00526619  0.9174892   0.39772583]
     [ 0.00000441 -0.39773137  0.91750191]]

You should find the matrix easy to work with if ``t`` is a single time,
but if ``t`` is a whole :ref:`array of times <date-arrays>`
then the resulting matrix will have a third dimension
with the same number of elements as the time vector.
NumPy provides no direct support
for rotation matrices with an extra dimension,
so avoid using NumPy’s multiplication operators.
Instead, use Skyfield utility functions:

.. testsetup::

    from numpy import array, identity
    R2 = identity(3)
    v = array([1,2,3])

.. testcode::

    from skyfield.functions import T, mxm, mxv

    T(R)        # reverse rotation matrix: frame → ICRS
    mxm(R, R2)  # matrix × matrix: combines rotations
    mxv(R, v)   # matrix × vector: rotates a vector
