
=============================
 Using Skyfield with AstroPy
=============================

.. currentmodule:: skyfield.api

The `AstroPy Project`_ is a sprawling community effort
to bring together a complete toolkit
for working astronomers and astrophysicists.
It can convert between dozens of different units,
wrap numeric vectors in its own array and table data types,
translate between a number of timescales and coordinate frames,
run image processing algorithms on astronomical images,
and query online sky catalogs.

Skyfield does not depend upon AstroPy,
but can represent its results using AstroPy units.
Here are the points of connection that Skyfield provides
between the two libraries:

1. You can provide an AstroPy ``Time`` value
   and Skyfield will translate it into its own representation.

   .. testcode::

      from astropy.time import Time
      from skyfield.api import load

      atime = Time('2010-01-01T00:00:00', scale='utc')
      print(atime)

      ts = load.timescale()
      t = ts.from_astropy(atime)
      print(t.utc_jpl())

   .. testoutput::

      2010-01-01T00:00:00.000
      A.D. 2010-Jan-01 00:00:00.0000 UTC

2. A skyfield
   :class:`~skyfield.positionlib.Barycentric`,
   :class:`~skyfield.positionlib.Astrometric`, or
   :class:`~skyfield.positionlib.Apparent`
   position and velocity
   can convert themselves into an AstroPy quantity
   using any linear and velocity units you specify.

   .. testcode::

      import astropy.units as u
      from skyfield.api import load

      planets = load('de421.bsp')
      earth = planets['earth']

      ts = load.timescale()
      t = ts.utc(1980, 1, 1)
      barycentric = earth.at(t)

      print(barycentric.position.to(u.au))
      print(barycentric.velocity.to(u.au / u.day))

   .. testoutput::

      [-0.16287311  0.88787399  0.38473904] AU
      [-0.01721258 -0.00279426 -0.0012121 ] AU / d

3. A Skyfield position can also return a complete AstroPy ``SkyCoord``
   object that couples the position vector with its reference frame.

   .. testcode::

      from astropy.coordinates import ICRS
      sc = barycentric.to_skycoord()
      print(sc)

   .. testoutput::

      <SkyCoord (ICRS): (x, y, z) in AU
          (-0.16287311, 0.88787399, 0.38473904)>

4. A Skyfield angle can express itself as an AstroPy quantity
   in any requested unit of angular measure.

   .. testcode::

      ra, dec, distance = barycentric.radec()
      declination = dec.to(u.deg)
      print('{0:0.03f}'.format(declination))

   .. testoutput::

      23.084 deg

.. _Astropy Project: http://docs.astropy.org/en/stable/

