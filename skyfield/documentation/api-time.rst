
======================
 API Reference — Time
======================

See the `time` guide for a careful and detailed discussion
of the several time scales used by astronomers,
and of how you can convert time to and from familiar time scales
like UTC and the worldwide time zones that are adapted from it.

Timescale, for building and converting times
============================================

.. currentmodule:: skyfield.timelib

.. autoclass:: Timescale
   :members:

The Time object
===============

.. autoclass:: Time
   :members:

   Four basic floating-point values
   can be directly accessed as attributes:

   .. attribute:: tai

      International Atomic Time (TAI) as a Julian date.

   .. attribute:: tt

      Terrestrial Time (TT) as a Julian date.

   .. attribute:: J

      Terrestrial Time (TT) as a floating point number of Julian years.

   .. attribute:: tdb

      Barycentric Dynamical Time (TDB) as a Julian date.

   .. attribute:: ut1

      Universal Time (UT1) as a Julian date.

   Two standard differences between time scales
   are also available as attributes:

   .. attribute:: delta_t

      The difference TT − UT1 measured in seconds.

   .. attribute:: dut1

      The difference UT1 − UTC measured in seconds.

   All of the other ways of expressing the time
   and converting it to typical human systems
   like UTC and world time zones
   are offered through methods:

.. _datetime: https://docs.python.org/3.5/library/datetime.html#datetime.datetime
.. _pytz: http://pythonhosted.org/pytz/
