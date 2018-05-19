
===============
 API Reference
===============

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
   Loader.path_to
   Loader.timescale
   Loader.tle

Time scales
===========

.. currentmodule:: skyfield.timelib

A Skyfield `Timescale` object is typically built
at the beginning of each program::

    from skyfield import api
    ts = api.load_timescale()

It downloads and parses the data tables necessary
to correctly convert between Universal Time
and the more stable time scales used by astronomers.

.. autosummary::

   Timescale.now
   Timescale.utc
   Timescale.tai
   Timescale.tt
   Timescale.tdb
   Timescale.from_astropy

Time objects
============

The `Time` class is Skyfield's way of representing
either a single time, or a whole array of times.
Each object has four floating point attributes
that present the time in several basic time scales.

========= ==================================================
``t.tai`` International Atomic Time (TAI) as a Julian date.
``t.tt``  Terrestrial Time (TT) as a Julian date.
``t.tdb`` Barycentric Dynamical Time (TDB) as a Julian date.
``t.ut1`` Universal Time (UT1) as a Julian date.
========= ==================================================

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

Vector Functions
================

.. currentmodule:: skyfield.vectorlib

The common API shared by planets, Earth locations, and Earth satellites.

.. autosummary::

   VectorFunction
   VectorFunction.__add__
   VectorFunction.__sub__
   VectorFunction.at

Planetary Ephemerides
=====================

.. currentmodule:: skyfield.jpllib

By downloading a `SpiceKernel` file,
Skyfield users can build vector functions
predicting the positions of the Moon, Sun, and planets.

.. autosummary::

   SpiceKernel
   SpiceKernel.comments
   SpiceKernel.names
   SpiceKernel.decode
   SpiceKernel.__getitem__

Topocentric Locations
=====================

.. currentmodule:: skyfield.toposlib

You can create a vector function
that computes the location of any position on the Earth’s surface.

.. autosummary::

   Topos

Earth Satellites
================

.. currentmodule:: skyfield.sgp4lib

By downloading TLE satellite element sets,
Skyfield users can build vector functions
that predict their positions.

.. autosummary::

   EarthSatellite

Stars and other distant objects
===============================

.. currentmodule:: skyfield.starlib

.. autosummary::
   :nosignatures:

   Star

Astronomical positions
======================

.. currentmodule:: skyfield.positionlib

The `ICRF` class serves as the base for all other positions classes,
which each share its axes but have more specific meanings.
Positions are usually returned by calling ``at(t)`` on a vector function
rather than being constructed manually.

.. autosummary::
   :nosignatures:

   ICRF
   Barycentric
   Astrometric
   Apparent
   Geocentric

Position methods and attributes
===============================

All position objects have three attributes
which provide access to their raw data:

================= ==================================
``ICRF.t``        The `Time` of the position.
``ICRF.position`` A `Distance` array giving x, y, z.
``ICRF.velocity`` A `Velocity` array giving ẋ, ẏ, ż.
================= ==================================

If a position lacks a velocity,
then the attribute is simply the value ``None``.

All positions support a basic set of methods:

.. autosummary::

   ICRF.distance
   ICRF.speed
   ICRF.radec
   ICRF.separation_from
   ICRF.ecliptic_position
   ICRF.ecliptic_latlon
   ICRF.galactic_position
   ICRF.galactic_latlon
   ICRF.from_altaz

Position methods specific to one class
======================================

.. autosummary::

   Barycentric.observe
   Astrometric.apparent
   Apparent.altaz
   Geocentric.subpoint

Units
=====

.. currentmodule:: skyfield.units

===================== ==================================================
`Distance`            Distance measure.
``Distance.au``       Astronomical Units.
``Distance.km``       Kilometers.
``Distance.m``        Meters.
`Velocity`            Velocity measure.
``Velocity.au_per_d`` Astronomical Units.
``Velocity.km_per_s`` Kilometers.
`Angle`               Angle measure.
``Angle.degrees``     Degrees of arc (360 in a complete circle).
``Angle.hours``       Hours of arc (24 in a complete circle).
``Angle.radians``     Radians (τ = 2π in a complete circle).
===================== ==================================================

All three kinds of quantity support one or more methods.

.. autosummary::

   Distance.to
   Velocity.to
   Angle.to
   Angle.hms
   Angle.signed_hms
   Angle.hstr
   Angle.dms
   Angle.signed_dms
   Angle.dstr
