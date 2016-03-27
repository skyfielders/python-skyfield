
==================================
 API Reference (work in progress)
==================================

.. currentmodule:: skyfield.jpllib

Ephemerides
===========

Skyfield users will usually download and open a `SpiceKernel` file
in a single step by calling `load()`.

.. autosummary::

   SpiceKernel
   SpiceKernel.decode
   SpiceKernel.__getitem__
   Body.at
   Body.geometry_of
   Body.topos
   Geometry.at

.. currentmodule:: skyfield.positionlib

Positions
=========

The `ICRF` class serves as the base for all other positions classes,
which each share its axes but have more specific meanings:

.. autosummary::
   :nosignatures:

   ICRF
   Barycentric
   Astrometric
   Apparent
   Geocentric

Position methods
================

.. autosummary::

   ICRF.distance
   ICRF.speed
   ICRF.radec
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
   Geocentric.observe

.. testsetup::

    from skyfield.positionlib import *
    from skyfield.api import load
    ts = load.timescale()

.. currentmodule:: skyfield.jpllib

.. testsetup::

   from skyfield import api
   de421 = api.load('de421.bsp')
   earth = de421['Earth']
   moon = de421['Moon']

.. autoclass:: SpiceKernel
   :members:

.. autoclass:: Body
   :members:

.. autoclass:: Geometry
   :members:

Generic ICRF position
=====================

.. currentmodule:: skyfield.positionlib

.. autoclass:: ICRF
   :members:

Position centered on the Solar System barycenter
================================================

.. autoclass:: Barycentric
   :members:

Astrometric position centered on an observer
============================================

.. autoclass:: Astrometric
   :members:

Apparent position centered on an observer
=========================================

.. autoclass:: Apparent
   :members:

.. autoclass:: Geocentric
   :members:

Timescale, for building and converting times
============================================

.. currentmodule:: skyfield.timelib

.. autoclass:: Timescale
   :members:

The Time object
===============

.. autoclass:: Time
   :members:
