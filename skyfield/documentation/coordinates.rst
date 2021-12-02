
=============
 Coordinates
=============

Add text to Positions:

at end of radec
If instead you are starting with ra/dec coordinates
and want to convert them to a Skyfield position,
see _.

And do that with other ra/dec? and altaz?

“Coordinates” are numbers that we use to name astronomical positions.
For any given position,
Skyfield can return at least a half-dozen
different coordinates,
each of which is simply a different but equivalent set of numbers
for designating the same point in the sky.

A few coordinate systems are so important
that they’ve already been explained in the :doc:`positions` chapter.
You can read about them there:

* :ref:`positions:Right Ascension and Declination: Astrometric`
* :ref:`positions:Right Ascension and Declination: Apparent`
* :ref:`positions:Azimuth and altitude from a geographic position`

This chapter introduces some further coordinate systems
that Skyfield supports.

Through the discussion that follows,
remember that coordinates are merely names for a position.
The choice of coordinate don’t change the position itself.
As explained in the :doc:`positions` chapter,
there’s a real physical difference between,
say, an astrometric position and an apparent position:
they are two distinct points in the sky.
But once you’ve used Skyfield to compute a particular position,
the choice of coordinate simply the choice
of what kind of label you want Skyfield to generate
for that exact point in the sky.

Cartesian coordinates versus Spherical coordinates
==================================================

A single coordinate system can be used to produce
two kinds of coordinate —
either Cartesian |xyz| coordinates,
or else two angles and a distance.
You will recall that both were illustrated
in the :doc:`positions` chapter.
There we generated both |xyz| coordinates and angular coordinates
that were measured against the axes
of the ICRS high-precision reference frame
that Skyfield uses internally:

.. testcode::

    # From the chapter on Positions: building
    # a position and asking for coordinates.

    from skyfield.api import load

    ts = load.timescale()
    t = ts.utc(2021, 12, 2, 14, 7)
    eph = load('de421.bsp')
    earth, mars = eph['earth'], eph['mars']
    astrometric = earth.at(t).observe(mars)

    # The position in Cartesian coordinates.

    x, y, z = astrometric.xyz.au
    print('Cartesian:')
    print('  x = {:.3f} au'.format(x))
    print('  y = {:.3f} au'.format(y))
    print('  z = {:.3f} au'.format(z))
    print()

    # The position in spherical coordinates.

    ra, dec, distance = astrometric.radec()
    print('Spherical:')
    print(' ', ra, 'right ascension')
    print(' ', dec, 'declination')
    print(' ', distance, 'distance')

.. testoutput::

    Cartesian:
      x = -1.521 au
      y = -1.800 au
      z = -0.772 au

    Spherical:
      15h 19m 15.08s right ascension
      -18deg 08' 37.0" declination
      2.47983 au distance

A note on terminology:
coordinates that measure simple linear distances
along three orthogonal axes
are called *Cartesian*, after their inventor René Descartes.
Whereas coordinates that measure two angles and a distance
are called *spherical*
because the two angles designate a position on the celestial sphere
just like the coordinates you might read off of the lines on a globe.
This makes spherical coordinates very convenient for observations.
Since a telescope doesn’t care whether an object is close or distant,
but just needs to know where to point,
all it needs are the two angles.

Each particular coordinate system, then, involves two steps.
First it must specify a reference frame,
the *x*-axis, *y*-axis, and *z*-axis
along which it measures its coordinates.
Second, it must specify a convention
for turning those coordinates into angles.
For example,
the coordinate system printed out by the above code
can be summarized as:

 | **Equatorial Coordinates**
 | *xy*-plane: Earth’s equator
 | *x*-axis: March equinox
 | *z*-axis: North celestial pole
 | ↕ Declination ±90° from equator toward poles
 | ↔ Right ascension 0\ |h|–24\ |h| east around the equator

The :doc:`positions` chapter gave you a second example of a coordinate system
when it introduced the concepts of altitude and azimuth.
Together they form a second coordinate system,
called *horizontal* because it is anchored to the horizon:

 | **Horizonal Coordinates**
 | *xy*-plane: Horizon
 | *x*-axis: North point along the horizon
 | *z*-axis: Zenith
 | ↕ Altitude ±90° above or below horizon
 | ↔ Azimuth 0°–360° measured clockwise from north
 | And is a LEFT HANDED coordinate system,
 | such that *y* is positive toward the east

.. testcode::

    from skyfield.api import load, wgs84
    bluffton = wgs84.latlon(+40.8939, -83.8917)
    astrometric = (earth + bluffton).at(t).observe(mars)
    alt, az, _ = astrometric.apparent().altaz()
    print(alt)  # 24deg, so partway up
    print(az)   # 147deg, so SE/SSE
    x, y, z = astrometric.frame_xyz(bluffton).au
    print('Cartesian:')
    print('  x = {:.3f} au'.format(x))  # W: negative middlest
    print('  y = {:.3f} au'.format(y))  # N: negative biggest
    print('  z = {:.3f} au'.format(z))  # positive but small  YES!

.. testoutput::

    24deg 24' 20.2"
    147deg 54' 27.9"
    Cartesian:
      x = -1.913 au
      y = 1.200 au
      z = 1.025 au

As we now turn to several other coordinate systems,
we will see this pattern repeated:
a reference frame that supports |xyz| coordinates,
plus a convention for expressing the vector’s direction
as a pair of angles.

Ecliptic and Galactic coordinates
==================================

Both of these coordinate systems
borrow the terms *latitude* and *longitude* from geography,
because they each measure angles
against a great circle and a north and south pole
that stand at right angles to that circle.

*Ecliptic coordinates* are measured from the plane of the Earth’s orbit.
They are useful
when making maps and diagrams of the Solar System
and when exploring the properties of orbits around the Sun,
because it places the orbits of the major planets
nearly flat against the *xy*-plane —
unlike right ascension and declination,
which twists the Solar System up at the 23° tilt of the Earth’s own axis.

 | **Ecliptic Coordinates**
 | *xy*-plane: Ecliptic plane
 | *x*-axis: March equinox
 | *z*-axis: North ecliptic pole
 | ↕ Latitude ±90° above or below the ecliptic
 | ↔ Longitude 0°–360° measured east from equinox

You might be tempted to ask
why we measure against the plane of the Earth’s orbit,
instead of averaging together all the planets
to compute the “invariable plane” of the whole Solar System
(to which the Earth’s orbit is inclined by something like 1.57°).
The answer is: precision.
We know the plane of the Earth’s orbit to many decimal places,
because the Earth carries all of our highest-precision observatories
along with it as it revolves.
Our estimate of the invariable plane, by contrast,
is a mere average
that changes every time we discover a new trans-Neptunian object
and — by small amounts — every time a new asteroid or comet is discovered.
So our Earth’s own orbit remains the best basis
for a coordinate system oriented to the Solar System.

To distinguish this latitude and longitude from the terrestrial ones,
it’s best to always call them *ecliptic latitude* and *ecliptic longitude*
with the word “ecliptic” always in front of them.
In the same way,
be specific and say *galactic latitude* and *galactic longitude*
when measuring angles relative to the plane and center of our galaxy —
or at least as best as we can approximate those
from our vantage point here deep inside the Orion Arm:

 | **Galactic Coordinates**
 | *xy*-plane: Galactic plane
 | *x*-axis: Galactic center
 | *z*-axis: North galactic pole
 | ↕ Latitude ±90° above galactic plane
 | ↔ Longitude 0°–360° east from galactic center

Astronomers have generated a series of more and more precise estimates
of our galaxy’s orientation over the past hundred years.
Skyfield uses the `IAU 1958 Galactic System II
<https://adsabs.harvard.edu/full/1960MNRAS.121..123B>`_,
which is believed to be accurate to within ±0.1°.

You can produce coordinates for either of these systems
by importing the ``ecliptic`` and ``galactic`` objects
from Skyfield’s frames library.
For example,
here’s how to compute ecliptic coordinates:

    from skyfield.framelib import ecliptic

    x, y, z = astrometric.frame_xyz(ecliptic)
    print('Cartesian:')
    print('  x =', x)
    print('  y =', y)
    print('  z =', z)

And here are the corresponding angles:

    lat, lon, distance = astrometric.frame_latlon(ecliptic)
    print('Spherical:')
    print(' ', lat, 'ecliptic latitude')
    print(' ', lon, 'ecliptic longitude')
    print(' ', distance, 'distance')

To produce galactic coordinates,
simply edit this code to use Skyfield’s `galactic` frame object instead.

///
But |xyz| coordinates are of limited use
when you want to point a telescope.
For that task you want spherical coordinates,
because they isolate the direction of the vector —
which is the only thing the telescope cares about —
from the vector’s length.

//all these can be applied to astrometric or apparent

Use the ``ecliptic``.

ECI and ECEF coordinates
========================

[Move over from positions.html]



.. TODO Make a full list of reference frames, including equatorial
   for the rare case where someone wants precessed/nuted xyz,
   and the planetary coordinate systems.

List of coordinate systems
==========================

horizontal

Turning coordinates into a position
===================================

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
