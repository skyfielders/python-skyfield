
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
   SpiceKernel.comments
   SpiceKernel.names
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

Position methods and attributes
===============================

.. autosummary::

   ICRF.t
   ICRF.position
   ICRF.velocity
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

    from pprint import pprint
    from skyfield.positionlib import *
    from skyfield.api import load
    ts = load.timescale()

.. currentmodule:: skyfield.jpllib

.. testsetup::

   from skyfield import api
   de421 = api.load('de421.bsp')
   earth = de421['Earth']
   moon = de421['Moon']
   mars = de421['Mars']

.. autoclass:: SpiceKernel
   :members:

.. autoclass:: Body
   :members:

.. autoclass:: Geometry()

   .. automethod:: at

Generic ICRF position
=====================

.. currentmodule:: skyfield.positionlib

.. autoclass:: ICRF
   :members:

   .. attribute:: t

      The `Time` coordinate of this position.

   .. attribute:: position

      The `Distance` coordinate as an (x, y, z) array.

   .. attribute:: velocity

      The `Velocity` coordinate as an (x, y, z) array.

      This attribute will have the value `None` if no velocity was
      specified for this position.

Position measured from the Solar System barycenter
==================================================

.. autoclass:: Barycentric
   :members:

Astrometric position relative to an observer
============================================

.. autoclass:: Astrometric
   :members:

Apparent position relative to an observer
=========================================

.. autoclass:: Apparent
   :members:

Geocentric position relative to the Earth
=========================================

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
