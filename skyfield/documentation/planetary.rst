
============================
 Planetary Reference Frames
============================

Skyfield has limited support for planetary reference frames
as defined in Jet Propulsion Lab data files.
It supports:

* Loading ``.tf`` and ``.tpc`` text files
  defining a wide array of planetary constants.

* Loading ``.bcp`` binary arrays
  that predict the orientation in space of Solar System bodies
  over a range of dates.

* Given a latitude and longitude on a Solar System body like the Moon that is
  (a) spherical and
  (b) whose orientation is defined by a data series in a ``.bcp`` file,
  computing that surface position’s location in space —
  allowing observations either from the point of view of that surface position,
  or observations of that surface position
  for other observers in the Solar System.

This leaves several features of such files still unsupported, though.
Skyfield:

* *Cannot* yet support a planetary reference frame
  whose orientation is described as a numeric series
  directly a text file.

* *Cannot* yet compute latitude and longitude positions
  on the surface of non-spherical objects,
  whose three axes in the text file are not the same length.

It is expected that support for these features will be added someday,
making Skyfield’s support for planetary constants complete.

Observing a Moon location
=========================

Here is how you would load up enough data
to predict where in the sky you would point a telescope
to see a particular latitude and longitude on the Moon —
in this example, the famous Aristarchus crater.

.. testcode::

    from skyfield.api import PlanetaryConstants, load

    ts = load.timescale()
    t = ts.utc(2019, 12, 20, 11, 5)

    eph = load('de421.bsp')
    earth, moon = eph['earth'], eph['moon']

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_text(load('pck00008.tpc'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')
    aristarchus = moon + pc.build_latlon_degrees(frame, 26.3, -46.8)

    apparent = earth.at(t).observe(aristarchus).apparent()
    ra, dec, distance = apparent.radec(epoch='date')
    print(ra)
    print(dec)

.. testoutput::

    13h 03m 22.96s
    -00deg 55' 27.3"

If your Moon location has a nonzero elevation
placing it above or below the Moon’s “sea level”,
you can provide ``build_latlon_degrees()``
with an extra ``elevation_m`` argument.

Observing from a Moon location
==============================

You can also ask Skyfield:
where in the sky would an astronaut standing on the Moon look
to see a particular Solar System body?
The answer can be provided in either right ascension and declination
coordinates against the background of stars,
or in altitude and azimuth measured against the astronaut’s horizon.

.. testcode::

    from skyfield.api import PlanetaryConstants, load

    ts = load.timescale()
    t = ts.utc(2019, 12, 20, 11, 5)

    eph = load('de421.bsp')
    earth, moon = eph['earth'], eph['moon']

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_text(load('pck00008.tpc'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')
    aristarchus = moon + pc.build_latlon_degrees(frame, 26.3, -46.8)

    apparent = aristarchus.at(t).observe(earth).apparent()
    ra, dec, distance = apparent.radec(epoch='date')
    print(ra)
    print(dec)

    alt, az, distance = apparent.altaz()
    print(alt, 'degrees above the horizon')
    print(az, 'degrees around the horizon from north')

.. testoutput::

    01h 03m 22.96s
    +00deg 55' 27.3"
    32deg 27' 09.7" degrees above the horizon
    118deg 12' 55.9" degrees around the horizon from north

Computing the sub-solar point on the Moon
=========================================

This works not only for the Sun,
but for any target body.
You can learn the Moon latitude and longitude directly beneath the body
by observing the target from the Moon’s center
and then asking the lunar reference frame
for the latitude and longitude.

.. testcode::

    sun = eph['Sun']

    p = moon.at(t).observe(sun).apparent()
    lat, lon, distance = p.frame_latlon(frame)
    lon_degrees = (lon.degrees + 180.0) % 360.0 - 180.0
    print('Sub-solar latitude: {:.1f} degrees'.format(lat.degrees))
    print('Sub-solar longitude: {:.1f} degrees'.format(lon_degrees))

.. testoutput::

    Sub-solar latitude: 0.3 degrees
    Sub-solar longitude: -104.9 degrees

Computing lunar libration
=========================

The Moon’s libration is expressed
as the latitude and longitude of the Moon location
that is currently nearest the Earth.
The convention seems to be that the simple geometric difference
between the Earth’s and Moon’s positions are used,
rather than the light-delayed position.
Thus:

.. testcode::

    p = (earth - moon).at(t)
    lat, lon, distance = p.frame_latlon(frame)
    lon_degrees = (lon.degrees + 180.0) % 360.0 - 180.0
    print('Libration in latitude: {:.3f}'.format(lat.degrees))
    print('Libration in longitude: {:.3f}'.format(lon_degrees))

.. testoutput::

    Libration in latitude: -6.749
    Libration in longitude: 1.520

The only subtlety is that the libration longitude
is not expressed as a number between 0° and 360°,
as would be more usual for longitude,
but instead as an offset positive or negative from zero,
which the above code accomplishes with some quick subtraction and modulo.

Computing a raw rotation matrix
===============================

If you are directly manipulating vectors,
you might simply want Skyfield to compute the NumPy rotation matrix
for rotating vectors from the ICRF into the frame of reference
of the Solar System body.
The ``frame`` object returned above
can return these matrices directly.
If given a single time ``t``,
the result will be a simple 3×3 matrix.

.. testcode::

    from skyfield.api import PlanetaryConstants, load

    ts = load.timescale()
    t = ts.utc(2019, 12, 20, 11, 5)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')

    R = frame.rotation_at(t)
    print(R.shape)

.. testoutput::

    (3, 3)

An array of times, by contrast,
will return an array of matrices
whose last dimension is as deep as the time vector is long.

.. testcode::

    t = ts.utc(2019, 12, 20, 11, range(5, 15))
    R = frame.rotation_at(t)
    print(t.shape)
    print(R.shape)

.. testoutput::

    (10,)
    (3, 3, 10)

The transpose ``R.T`` of the rotation matrix
can be used to rotate vectors
that are already in the reference frame of the body
back into a standard ICRF vector.
