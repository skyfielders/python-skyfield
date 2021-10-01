
======================================
 API Reference — Geographic Locations
======================================

.. currentmodule:: skyfield.toposlib

.. autodata:: wgs84
.. autodata:: iers2010

.. autoclass:: Geoid
   :members:

.. autoclass:: GeographicPosition
   :members:

   .. attribute:: model

      The :class:`Geoid`, like WGS84 or IERS2010,
      that this position uses to map longitude, latitude, and elevation
      to a three-dimensional Cartesian position.

   .. attribute:: latitude

      An :class:`~skyfield.units.Angle` object
      specifying the latitude of the topocentric position.
      The north pole has latitude +90°.

   .. attribute:: longitude

      An :class:`~skyfield.units.Angle` object
      specifying the longitude of the topocentric position.
      East is positive.

   .. attribute:: elevation

      A :class:`~skyfield.units.Distance` object
      specifying an elevation above (positive) or below (negative)
      the mean sea level of the Earth ellipsoid
      specified by this position’s :attr:`~GeographicPosition.model`.

   .. attribute:: itrs_xyz

      A :class:`~skyfield.units.Distance` object
      giving the spatial |xyz| coordinates of this location
      in the ITRS Earth-centered Earth-fixed (“ECEF”) reference frame.

   .. attribute:: center

      The integer 399,
      which identifies this topocentric position’s vector
      as having its origin at the center of the Earth.

   .. method:: at(t)

      Return the position of this Earth location at time ``t``.
