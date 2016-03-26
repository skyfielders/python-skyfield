
==================================
 API Reference (work in progress)
==================================

.. currentmodule:: skyfield.positionlib

.. autosummary::
   :nosignatures:

   ICRF
   Barycentric
   Astrometric
   Apparent

.. autosummary::

   ICRF.distance
   ICRF.radec
   ICRF.radec
   ICRF.ecliptic_position
   ICRF.ecliptic_latlon
   ICRF.galactic_position
   ICRF.galactic_latlon

.. testsetup::

    from skyfield.positionlib import *

Generic ICRF position
=====================

.. autoclass:: skyfield.positionlib.ICRF

   .. automethod:: skyfield.positionlib.ICRF.distance()

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
