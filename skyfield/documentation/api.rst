
==================================
 API Reference (work in progress)
==================================

.. currentmodule:: skyfield.jpllib

Ephemerides
===========

You generally

.. autosummary::

   SpiceKernel
   SpiceKernel.decode
   SpiceKernel.__getitem__
   Body.at
   Body.topos

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

Position methods, general
=========================

.. autosummary::

   ICRF.distance
   ICRF.speed
   ICRF.radec
   ICRF.ecliptic_position
   ICRF.ecliptic_latlon
   ICRF.galactic_position
   ICRF.galactic_latlon
   ICRF.from_altaz

Position methods, specific
==========================

.. autosummary::

   Barycentric.observe
   Astrometric.apparent
   Apparent.altaz
   Geocentric.observe

.. testsetup::

    from skyfield.positionlib import *
    from skyfield.api import load
    ts = load.timescale()

.. autoclass:: skyfield.jpllib.SpiceKernel

.. autoclass:: skyfield.jpllib.Body

Generic ICRF position
=====================

.. autoclass:: skyfield.positionlib.ICRF
   :members:

Position centered on the Solar System barycenter
================================================

.. autoclass:: skyfield.positionlib.Barycentric
   :members:

Astrometric position centered on an observer
============================================

.. autoclass:: skyfield.positionlib.Astrometric
   :members:

Apparent position centered on an observer
=========================================

.. autoclass:: skyfield.positionlib.Apparent
   :members:

.. autoclass:: skyfield.positionlib.Geocentric
   :members:

Timescale, for building and converting times
============================================

.. autoclass:: skyfield.api.Timescale
   :members:

The Time object
===============

.. autoclass:: skyfield.api.Time
   :members:
