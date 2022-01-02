
===============
 API Reference
===============

Quick links to the sections below:

.. contents:: :local:

Version
=======

Skyfield offers a tuple ``skyfield.VERSION``
that lets your code determine the installed version of Skyfield.

::

    import skyfield
    print(skyfield.VERSION)

See :ref:`Skyfield VERSION constant`.

Opening files
=============

.. currentmodule:: skyfield.iokit

::

    # File you already have.

    from skyfield.api import load_file
    planets = load_file('~/Downloads/de405.bsp')

.. autosummary::

   load_file

::

    # File you want Skyfield to download automatically.

    from skyfield.api import load
    ts = load.timescale()
    planets = load('de405.bsp')

.. autosummary::

   Loader
   Loader.build_url
   Loader.days_old
   Loader.download
   Loader.path_to
   Loader.timescale
   Loader.tle_file

.. _api-Timescale:

.. _timescale-summary:

Time scales
===========

.. currentmodule:: skyfield.timelib

A script will typically start by building a single Skyfield `Timescale`
to use for all date and time conversions:

.. testcode::

    from skyfield import api
    ts = api.load.timescale()

Its methods are:

.. autosummary::

   Timescale.now
   Timescale.from_datetime
   Timescale.from_datetimes
   Timescale.utc
   Timescale.tai
   Timescale.tai_jd
   Timescale.tt
   Timescale.tt_jd
   Timescale.J
   Timescale.tdb
   Timescale.tdb_jd
   Timescale.ut1
   Timescale.ut1_jd
   Timescale.from_astropy
   Timescale.linspace

.. _api-Time:

Time objects
============

The `Time` class is Skyfield's way of representing
either a single time, or a whole array of times.
The same time can be represented in several different time scales.

========= =====================================================
``t.tai`` International Atomic Time (TAI) as a Julian date.
``t.tt``  Terrestrial Time (TT) as a Julian date.
``t.J``   Terrestrial Time (TT) as floating point Julian years.
``t.tdb`` Barycentric Dynamical Time (TDB) as a Julian date.
``t.ut1`` Universal Time (UT1) as a Julian date.
========= =====================================================

A couple of offsets between time scales are also available.

============= ================================
``t.delta_t`` Difference TT − UT1 in seconds.
``t.dut1``    Difference UT1 − UTC in seconds.
============= ================================

Other time scales and conversions are available through its methods.

.. autosummary::

   Time.utc_jpl
   Time.utc_iso
   Time.utc_strftime
   Time.utc_datetime
   Time.utc_datetime_and_leap_second
   Time.astimezone
   Time.astimezone_and_leap_second
   Time.toordinal
   Time.tai_calendar
   Time.tt_calendar
   Time.tdb_calendar
   Time.ut1_calendar
   Time.tai_strftime
   Time.tt_strftime
   Time.tdb_strftime
   Time.ut1_strftime
   Time.M
   Time.MT
   Time.gmst
   Time.gast
   Time.nutation_matrix
   Time.precession_matrix

Time utilities
==============

.. autosummary::

   compute_calendar_date

Vector functions
================

.. currentmodule:: skyfield.vectorlib

The common API shared by planets, Earth locations, and Earth satellites.

.. autosummary::

   VectorFunction
   VectorFunction.at

Either adding two vector functions ``v1 + v2`` or subtracting them ``v1 - v2``
produces a new function of time that, when invoked with ``.at(t)``,
returns the sum or difference of the vectors returned by the two functions.

Planetary ephemerides
=====================

.. currentmodule:: skyfield.jpllib

By downloading a `SpiceKernel` file,
Skyfield users can build vector functions
predicting the positions of the Moon, Sun, and planets.
See :doc:`planets`.

.. autosummary::

   SpiceKernel
   SpiceKernel.close
   SpiceKernel.comments
   SpiceKernel.names
   SpiceKernel.decode

Kernels also support lookup using the Python ``kernel['Mars']`` syntax,
in which case they return a function of time
that returns vectors from the Solar System barycenter to the named body.

Planetary magnitudes
====================

.. autofunction:: skyfield.magnitudelib.planetary_magnitude

Planetary reference frames
==========================

.. currentmodule:: skyfield.planetarylib

.. autosummary::
   :nosignatures:

   PlanetaryConstants
   Frame

Almanac
=======

.. currentmodule:: skyfield.almanac

Routines to search for events like sunrise, sunset, and Moon phase.

.. autosummary::

   seasons
   moon_phase
   moon_phases
   moon_nodes
   oppositions_conjunctions
   meridian_transits
   sunrise_sunset
   dark_twilight_day
   risings_and_settings

.. currentmodule:: skyfield.eclipselib

.. autosummary::

   lunar_eclipses

Geographic locations
====================

.. currentmodule:: skyfield.toposlib

.. _topocentric-locations:

Skyfield supports two Earth datums
for translating between latitude/longitude and Cartesian coordinates.
They each use a slightly different estimate of the Earth’s oblateness.
The most popular is WGS84, which is used by the world’s GPS devices:

* `wgs84`
* `iers2010`

Each datum offers a method for taking a latitude and longitude
and returning a `GeographicPosition` that knows its position in space:

.. autosummary::

   Geoid.latlon

Going in the other direction,
there are several methods for converting an existing Skyfield position
into latitude, longitude, and height:

.. autosummary::

   Geoid.latlon_of
   Geoid.height_of
   Geoid.geographic_position_of
   Geoid.subpoint_of

Once you have used either of the above approaches
to build a `GeographicPosition`,
it offers several methods:

.. autosummary::

   GeographicPosition.at
   GeographicPosition.lst_hours_at
   GeographicPosition.refract
   GeographicPosition.rotation_at

Kepler orbits
=============

See :doc:`kepler-orbits`
for computing the positions of comets, asteroids, and other minor planets.

Earth satellites
================

.. currentmodule:: skyfield.sgp4lib

By downloading TLE satellite element sets,
Skyfield users can build vector functions
that predict their positions.
See :doc:`earth-satellites`.

.. autosummary::

   EarthSatellite
   EarthSatellite.from_satrec
   TEME

Stars and other distant objects
===============================

.. currentmodule:: skyfield.starlib

.. autosummary::
   :nosignatures:

   Star

Astronomical positions
======================

.. currentmodule:: skyfield.positionlib

The `ICRF` three-dimensional position vector serves as the base class
for all of the following position classes.  Each class represents an
|xyz| ``.position`` and ``.velocity`` vector oriented to the axes of
the International Celestial Reference System (ICRS),
an inertial system that’s an update to J2000
and that does not rotate with respect to the universe.

.. autosummary::
   :nosignatures:

   ICRF
   Barycentric
   Astrometric
   Apparent
   Geocentric

Positions are usually generated by the ``at(t)`` method of a vector
function, rather than being constructed manually.  But you can also
build a position directly from a raw vector, or from right ascension and
declination coordinates with
:func:`~skyfield.positionlib.position_of_radec()`.

.. autosummary::

   position_of_radec

All position objects offer five basic attributes:

.. PAT START

============= ========================================
``.position`` An |xyz| `Distance`.
``.velocity`` An |xyz| `Velocity`, or ``None``.
``.t``        The `Time` of the position, or ``None``.
``.center``   Body the vector is measured from.
``.target``   Body the vector is measured to.
============= ========================================

.. PAT END

All positions support these methods:

.. autosummary::

   ICRF.distance
   ICRF.speed
   ICRF.radec
   ICRF.hadec
   ICRF.separation_from
   ICRF.frame_xyz
   ICRF.frame_xyz_and_velocity
   ICRF.frame_latlon
   ICRF.frame_latlon_and_rates
   ICRF.from_time_and_frame_vectors
   ICRF.phase_angle
   ICRF.fraction_illuminated
   ICRF.is_sunlit
   ICRF.from_altaz

In addition to the methods above,
several subclasses of the base position class
provide unique methods of their own:

.. autosummary::

   Barycentric.observe
   Astrometric.apparent
   Apparent.altaz

Reference frames
================

.. autosummary::

   skyfield.framelib.true_equator_and_equinox_of_date
   skyfield.framelib.itrs
   skyfield.framelib.ecliptic_frame
   skyfield.framelib.ecliptic_J2000_frame
   skyfield.framelib.galactic_frame

Constellations
==============

.. autofunction:: skyfield.api.load_constellation_map
.. autofunction:: skyfield.api.load_constellation_names
.. autofunction:: skyfield.data.stellarium.parse_constellations
.. autofunction:: skyfield.data.stellarium.parse_star_names

Searching
=========

.. currentmodule:: skyfield.searchlib

.. autofunction:: find_discrete()
.. autofunction:: find_maxima()
.. autofunction:: find_minima()

Osculating orbital elements
===========================

This routine returns osculating orbital elements for an object’s
instantaneous position and velocity.

.. currentmodule:: skyfield.elementslib

.. autosummary::

   osculating_elements_of

================================================== ============================
``OsculatingElements.apoapsis_distance``           Distance object
``OsculatingElements.argument_of_latitude``        Angle object
``OsculatingElements.argument_of_periapsis``       Angle object
``OsculatingElements.eccentric_anomaly``           Angle object
``OsculatingElements.eccentricity``                numpy.ndarray
``OsculatingElements.inclination``                 Angle object
``OsculatingElements.longitude_of_ascending_node`` Angle object
``OsculatingElements.longitude_of_periapsis``      Angle object
``OsculatingElements.mean_anomaly``                Angle object
``OsculatingElements.mean_longitude``              Angle object
``OsculatingElements.mean_motion_per_day``         Angle object
``OsculatingElements.periapsis_distance``          Distance object
``OsculatingElements.periapsis_time``              Time object
``OsculatingElements.period_in_days``              numpy.ndarray
``OsculatingElements.semi_latus_rectum``           Distance object
``OsculatingElements.semi_major_axis``             Distance object
``OsculatingElements.semi_minor_axis``             Distance object
``OsculatingElements.time``                        Time object
``OsculatingElements.true_anomaly``                Angle object
``OsculatingElements.true_longitude``              Angle object
================================================== ============================

Units
=====

.. currentmodule:: skyfield.units

===================== ==================================================
`Distance`            Distance
`Velocity`            Velocity
`Angle`               Angle
`AngleRate`           Rate at which an angle is changing
===================== ==================================================

All three kinds of quantity support one or more methods.

.. autosummary::

   Distance.au
   Distance.km
   Distance.m
   Distance.length
   Distance.light_seconds
   Distance.to
   Velocity.au_per_d
   Velocity.km_per_s
   Velocity.m_per_s
   Velocity.to
   Angle.radians
   Angle.hours
   Angle.degrees
   Angle.arcminutes
   Angle.arcseconds
   Angle.mas
   Angle.to
   Angle.hms
   Angle.signed_hms
   Angle.hstr
   Angle.dms
   Angle.signed_dms
   Angle.dstr
   AngleRate.radians
   AngleRate.degrees
   AngleRate.arcminutes
   AngleRate.arcseconds
   AngleRate.mas
   Rate.per_day
   Rate.per_hour
   Rate.per_minute
   Rate.per_second

Trigonometry
============

.. currentmodule:: skyfield.trigonometry

.. autosummary::

   position_angle_of
