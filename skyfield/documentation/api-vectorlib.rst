
============================================================
 API Reference — Generating Positions From Vector Functions
============================================================


.. testsetup::

   from pprint import pprint
   from skyfield.api import load
   de421 = load('de421.bsp')
   earth = de421['Earth']
   moon = de421['Moon']

.. currentmodule:: skyfield.jpllib

JPL .bsp ephemeris files
========================

.. autoclass:: SpiceKernel
   :members:
