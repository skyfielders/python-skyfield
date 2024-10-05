
======================================
 API Reference — Geographic Locations
======================================

.. currentmodule:: skyfield.toposlib

.. autodata:: wgs84
.. autodata:: iers2010

.. autoclass:: Geoid
   :members:
   :exclude-members: subpoint

   .. method:: subpoint(…)

      .. deprecated:: 1.40

         Renamed to `geographic_position_of()`.

.. autoclass:: GeographicPosition
   :members:

   .. attribute:: model

      The :class:`Geoid`, like WGS84 or IERS2010,
      that this position uses to map latitude, longitude, and elevation
      to a three-dimensional Cartesian position.

   .. attribute:: latitude

      An :class:`~skyfield.units.Angle` specifying latitude;
      the north pole has latitude +90°.

   .. attribute:: longitude

      An :class:`~skyfield.units.Angle` specifying longitude;
      east is positive, west is negative.

   .. attribute:: elevation

      A :class:`~skyfield.units.Distance` specifying elevation
      above (positive) or below (negative)
      the surface of the Earth ellipsoid
      specified by this position’s :attr:`~GeographicPosition.model`.

   .. attribute:: itrs_xyz

      A :class:`~skyfield.units.Distance` object
      giving the spatial |xyz| coordinates of this location
      in the ITRS Earth-centered Earth-fixed (“ECEF”) reference frame.

   .. attribute:: center

      The integer 399,
      which identifies this position as geocentric:
      its |xyz| coordinates are measured from the Earth’s center.

   .. method:: at(t)

      Return the position of this Earth location at time ``t``.

.. autoclass:: ITRSPosition
   :members:

   This |xyz| vector has no knowledge
   of standard geoids, latitude, or longitude,
   but is convenient if you already know
   the rectangular coordinates of a target’s location:

   .. testcode::

    from skyfield.api import Distance
    from skyfield.toposlib import ITRSPosition

    d = Distance(km=[-3918, -1887, 5209])
    p = ITRSPosition(d)

   .. method:: at(t)

      Return the GCRS position of this ITRS coordinate at time ``t``.
