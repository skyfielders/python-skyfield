=======================
 API Reference — Units
=======================

When you ask positions to return distances, velocities, or angles,
they return instances of the following classes.

.. currentmodule:: skyfield.units

.. autoclass:: Distance
   :members:

   .. attribute:: au

      The distance in astronomical units.

   .. attribute:: km

      The distance in kilometers.

   .. attribute:: m

      The distance in meters.

.. autoclass:: Velocity
   :members:

   .. attribute:: au_per_d

      The velocity in astronomical units per day.

   .. attribute:: km_per_s

      The velocity in kilometers per second.

.. autoclass:: Angle
   :members:

   .. attribute:: degrees

      The angle in degrees (360° in a circle).

   .. attribute:: hours

      The angle in hours (24 hours in a circle).

   .. attribute:: radians

      The angle in radians (𝜏 = 2𝜋 in a circle)
