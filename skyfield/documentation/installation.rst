
==============
 Installation
==============

.. currentmodule:: skyfield.api

Skyfield has only a single binary dependency,
the `NumPy <http://www.numpy.org/>`_ vector library,
and is designed to install cleanly with a single invocation
of the standard Python package tool::

    pip install skyfield

If trying to install Skyfield gives you errors about NumPy,
there are several other ways to get NumPy installed:

* It is best to simply use
  `a science-ready version of Python
  <http://www.scipy.org/install.html#scientific-python-distributions>`_
  that comes with NumPy built-in,
  like the `Anaconda <http://continuum.io/downloads>`_
  distribution.

* | There are several approaches described in the `SciPy install instructions <http://www.scipy.org/install.html>`_.

* You can download and run an official `NumPy installer
  <https://sourceforge.net/projects/numpy/files/NumPy/>`_.

* You can try to get the plain ``pip``-powered install working
  by giving your system a functioning C compiler
  that matches the compiler used to build Python.
  On Windows with Python 2.7, try installing the free
  `Visual Studio Express 2008 <http://go.microsoft.com/?linkid=7729279>`_.
  Mac users should install the “Xcode Command Line Tools”
  to give ``pip`` the superpower of being able to build and install
  binary Python dependencies like NumPy.

Note that Skyfield has not reached version 1.0 yet,
so tweaks to the API are still possible.
Read the :ref:`changelog` below to learn about recent adjustments.
You can protect your project from any abrupt API changes
by pinning a specific version of Skyfield
in your ``requirements.txt`` or ``setup.py`` or install instructions::

    skyfield==0.5

By preventing Skyfield from being upgraded
until you are ready to advance the version number yourself,
you can avoid being blindsided by improvements that take place
between now and the eventual 1.0 version.
If you find any problems or would like to suggest an improvement,
simply create an issue on the project’s GitHub page:

    https://github.com/brandon-rhodes/python-skyfield

Good luck!

.. _changelog:

Change Log
==========

0.6
---

* Skyfield now automatically generates an estimate for `delta_t` if the
  user does not supply their own `delta_t` when specifying a date.  This
  should make altitude and azimuth measurements much more precise.

* The leap-second table has been updated to include 2015 July 1.

* Both ecliptic and galactic coordinates are now supported.

0.5
---

* Skyfield has dropped the 16-megabyte JPL ephemeris DE421 as an install
  dependency, since users might choose another ephemeris, or might not
  need one at all.  You now ask for a SPICE ephemeris to be downloaded
  at runtime with a call like ``planets = load('de421.bsp')``.

* Planets are no longer offered as magic attributes, but are looked up
  through the square bracket operator.  So instead of typing
  ``planets.mars`` you should now type ``planets['mars']``.  You can run
  ``print(planets)`` to learn which bodies an ephemeris supports.

* | Ask for planet positions with ``body.at(jd)`` instead of ``body(jd)``.

* Per IAU 2012 Resolution B2, Skyfield now uses lowercase *au* for the
  astronomical unit, and defines it as exactly 149 597 870 700 meters.
  While this API change is awkward for existing users, I wanted to make
  the change while Skyfield is still pre-1.0.  If this breaks a program
  that you already have running, please remember that a quick ``pip``
  ``install`` ``skyfield==0.4`` will get you up and running again until
  you have time to edit your code and turn ``AU`` into ``au``.

0.4
---

* To prevent confusion, the Julian date :meth:`~JulianDate.astimezone()`
  and :meth:`~JulianDate.utc_datetime()` methods
  have been changed to return only a ``datetime`` object.
  If you also need a leap second flag returned,
  call the new methods :meth:`~JulianDate.astimezone_and_leap_second()`
  and :meth:`~JulianDate.utc_datetime_and_leap_second()`.

0.3
---

* The floating-point values of an angle
  ``a.radians``, ``a.degrees``, and ``a.hours``
  are now attributes instead of method calls.
