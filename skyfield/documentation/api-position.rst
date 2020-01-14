
========================================
 API Reference — Astronomical Positions
========================================

See `positions` for a detailed guide
to these various kinds of position that Skyfield can compute,
and to the selection of coordinate systems
that can be used to express them.

.. testsetup::

   from __future__ import print_function
   from skyfield.api import load
   from skyfield.positionlib import ICRF
   ts = load.timescale()
   de421 = load('de421.bsp')
   earth = de421['Earth']
   mars = de421['Mars']

.. currentmodule:: skyfield.positionlib

Generic ICRF position
=====================

.. autoclass:: ICRF
   :members:
   :exclude-members: ecliptic_position, galactic_position

   .. attribute:: t

      The `Time` coordinate of this position.

   .. attribute:: position

      The `Distance` coordinate as an (x, y, z) array.

   .. attribute:: velocity

      The `Velocity` coordinate as an (x, y, z) array.

      This attribute will have the value ``None`` if no velocity was
      specified for this position.

Position measured from the Solar System barycenter
==================================================

.. autoclass:: Barycentric
   :members:

   This class inherits the methods of is parent class :class:`ICRF` as
   well as the orientation of its axes in space.

Astrometric position relative to an observer
============================================

.. autoclass:: Astrometric
   :members:

   This class inherits the methods of is parent class :class:`ICRF` as
   well as the orientation of its axes in space.

Apparent position relative to an observer
=========================================

.. autoclass:: Apparent
   :members:

   This class inherits the methods of is parent class :class:`ICRF` as
   well as the orientation of its axes in space.

Geocentric position relative to the Earth
=========================================

.. autoclass:: Geocentric
   :members:

   This class inherits the methods of is parent class :class:`ICRF` as
   well as the orientation of its axes in space.

Geometric instantaneous position between two objects
====================================================

.. autoclass:: Geometric
   :members:

   This class inherits the methods of is parent class :class:`ICRF` as
   well as the orientation of its axes in space.

Building a position from right ascension and declination
========================================================

.. autofunction:: position_from_radec
