
=======================================
 API Reference — Topocentric Locations
=======================================

.. currentmodule:: skyfield.toposlib

.. autoclass:: Topos
   :members:

   .. attribute:: latitude

      An :class:`~skyfield.units.Angle` object
      specifying the latitude of the topocentric position.

   .. attribute:: longitude

      An :class:`~skyfield.units.Angle` object
      specifying the longitude of the topocentric position.

   .. attribute:: elevation

      A :class:`~skyfield.units.Distance` object
      specifying the elevation of the topocentric position
      above mean sea level on a WGS-84 globe.

   .. attribute:: center

      The integer 399,
      which identifies this topocentric position’s vector
      as having its origin at the center of the Earth.
