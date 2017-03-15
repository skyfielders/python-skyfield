
==================================
 API Reference — Vector Functions
==================================

.. currentmodule:: skyfield.vectorlib

The Skyfield API is build atop *vector functions*
that take time as input and produce a position vector.
You can create vector functions for the Earth, Moon, Sun, planets,
and Earth satellites,
and can combine them using addition and subtraction.

.. autoclass:: VectorFunction
   :members:

   .. attribute:: center

      The Solar System object from which this vector is measured.
      Often this is an integer code like ``399`` for the Earth,
      ``3`` for the center of gravity of the Earth-Moon system,
      or ``0`` for the very center of the Solar System itself,
      though it might also be a specific object
      like a :class:`~skyfield.toposlib.Topos` on the Earth’s surface
      or an :class:`~skyfield.sgp4lib.EarthSatellite` in orbit around it.

   .. attribute:: target

      Using the same set of possible values as the ``center``,
      this attribute names the target to which the vector is pointing.
      The vector, then, is the three-dimensional difference
      between the position of the center and that of the target.

   .. attribute:: vf1 + vf2

      Return a new vector function whose ``at(t)``, when called,
      computes the sum of the original vectors ``vf1`` and ``vf2``.
      This will raise an error
      unless the ``target`` where one of the two vectors ends
      is the same as the ``center`` from which the other vector starts.

   .. attribute:: vf1 - vf2

      Return a new vector function whose ``at(t)``, when called,
      computes where the ``target`` of ``vf1`` will be positioned
      relative to the ``target`` of the subtracted ``vf2``.
      Note that this will be an instantaneous vector,
      uncorrected for the amount of time light takes
      to travel from one target to the other.
      This raises an error
      unless the two vectors share the same ``center``.
