
=======================================
 API Reference — Topocentric Locations
=======================================

.. currentmodule:: skyfield.toposlib

.. autodata:: wgs84
.. autodata:: iers2010

.. autoclass:: Geoid
   :members:

.. autoclass:: GeographicPosition
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

   .. attribute:: itrs_position

      A :class:`~skyfield.units.Distance` object
      giving the spatial x,y,z coordinates of this location
      in the ITRS Earth-centered Earth-fixed (“ECEF”) reference frame.

   .. attribute:: center

      The integer 399,
      which identifies this topocentric position’s vector
      as having its origin at the center of the Earth.

   .. method:: at(t)

      Return the position of this Earth location at time ``t``.
