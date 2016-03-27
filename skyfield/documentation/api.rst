
==================================
 API Reference (work in progress)
==================================

.. currentmodule:: skyfield.positionlib

Positions
=========

There is a single base class `ICRF` for positions,
along with several more specific classes
that inherit its methods:

.. autosummary::
   :nosignatures:

   ICRF
   Barycentric
   Astrometric
   Apparent
   Geocentric

The base methods shared by all position objects are:

.. autosummary::

   ICRF.distance
   ICRF.speed
   ICRF.radec
   ICRF.ecliptic_position
   ICRF.ecliptic_latlon
   ICRF.galactic_position
   ICRF.galactic_latlon
   ICRF.from_altaz

.. autosummary::

   

.. testsetup::

    from skyfield.positionlib import *
    from skyfield.api import load
    ts = load.timescale()

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
