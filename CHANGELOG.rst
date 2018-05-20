Changelog
=========

.. currentmodule:: skyfield.positionlib

1.4 — 2018 May 20
-----------------

* You can now specify the distance to an object when generating a
  position from altitude and azimuth coordinates.
  `#158 <https://github.com/skyfielders/python-skyfield/issues/158>`_

* The dictionary of satellites returned when you read a TLE file
  now supports lookup by integer satellite ID, not just by name,
  and now knows how to parse TLE files from Space-Track.
  `#163 <https://github.com/skyfielders/python-skyfield/issues/163>`_
  `#167 <https://github.com/skyfielders/python-skyfield/issues/167>`_

* Star coordinates can now be offered for any epoch, not just J2000.
  `#166 <https://github.com/skyfielders/python-skyfield/issues/166>`_

1.3 — 2018 April 15
-------------------

* Geocentric coordinates now have a
  :meth:`~skyfield.positionlib.Geocentric.subpoint()`
  method that computes the latitude and longitude
  of the point beneath that body.

* All of the ``Timescale`` time constructor methods now accept arrays.

* Emergency fix to stop Skyfield
  from endlessly downloading new copies of ``deltat.preds``,
  since the file has gone out of date at the USNO site.

* Fixed ability of a :class:`~skyfield.starlib.Star`
  to be initialized with a tuple that breaks units into minutes and seconds
  (broke in version 1.2).

* Issues fixed:
  `#170 <https://github.com/skyfielders/python-skyfield/issues/170>`_
  `#172 <https://github.com/skyfielders/python-skyfield/issues/172>`_

1.2 — 2018 March 29
-------------------

* The documentation now describes
  how to create an excerpt of a large JPL ephemeris
  without downloading the entire file.
  Several Skyfield tests now run much faster
  because they use an ephemeris excerpt instead of waiting for a download.

* For ``load_file()`` a leading ``~`` now means “your home directory”.

* You can now initialize a velocity from kilometers per second
  with ``Velocity(km_per_s=...)``.

* Empty time and angle objects no longer raise an exception when printed.
  (Thanks, JoshPaterson!)

* Issues fixed:
  `#160 <https://github.com/skyfielders/python-skyfield/issues/160>`_
  `#161 <https://github.com/skyfielders/python-skyfield/issues/161>`_
  `#162 <https://github.com/skyfielders/python-skyfield/issues/162>`_

1.1 — 2018 January 14
---------------------

* Positions can now be converted to AstroPy with
  :meth:`~skyfield.positionlib.ICRF.to_skycoord()`.

* You can now provide a timescale of your own to an
  :meth:`~skyfield.sgp4lib.EarthSatellite`
  instead of having it trying to load one itself.

* Downloaded files are no longer marked as executable on Windows.

* A friendly error message, rather than an obscure traceback, is now
  returned if you try converting a position to alt/az coordinates but
  the position was not measured from a position on the Earth’s surface.

1.0 — 2017 March 15
-------------------

* Brought the core API to maturity: replaced the narrow concept of
  building a “body” from several ephemeris segments with the general
  concept of a vector function that is the sum of several simpler vector
  functions.

* Added support for adding and subtracting vector functions.

* Deprecated the Earth ``topos()`` method in favor of vector addition.

* Deprecated the Earth ``satellite()`` method in favor of vector addition.

* Deprecated the body ``geometry_of()`` method in favor of vector subtraction.

* Celestrak satellite files can now be opened with ``load.tle(url_or_filename)``.

0.9.1 — 2016 December 10
------------------------

* Attempted to speed up Earth satellite calculations by caching a single
  time scale object instead of creating a new one each time.

* Fixed a possible divide-by-zero error when applying deflection to an
  apparent position.

0.9
---

* The ``observe()`` method of an observer on the Earth’s surface now
  correctly accounts for the way that the Earth’s gravity will deflect
  the apparent position of objects that are not exactly overhead,
  bringing Skyfield’s agreement with the Naval Observatory’s NOVAS
  library to within half a milliarcsecond.

* The time method ``tt_calendar()`` method no longer raises a
  ``TypeError`` when its value is an array.

* Running ``repr()`` on a ``Time`` array now produces a more compact
  string that only mentions the start and end of the time period.

* The ``api.load()`` call no longer attempts to animate a progress bar
  if the user is running it under IDLE, which would try to accumulate
  the updates as a single long line that eventually hangs the window.

0.8
---

* Added an `api` document to the project, in reverent imitation of the
  `Pandas API Reference`_ that I keep open in a browser tab every time I
  am using the Pandas library.

* New method `ICRF.separation_from()` computes the angular separation
  between two positions.

* Fixed ``==`` between `Time` objects and other unrelated objects so
  that it no longer raises an exception.

0.7
---

* Introduced the ``Timescale`` object with methods ``utc()``, ``tai()``,
  ``tt()``, and ``tdb()`` for building time objects, along with a
  ``load.timescale()`` method for building a new ``Timescale``.  The
  load method downloads ΔT and leap second data from official data
  sources and makes sure the files are kept up to date.  This replaces
  all former techniques for building and specifying dates and times.

* Renamed ``JulianDate`` to ``Time`` and switched from ``jd`` to ``t``
  as the typical variable used for time in the documentation.

* Deprecated timescale keyword arguments like ``utc=(…)`` for both the
  ``Time`` constructor and also for all methods that take time as
  an argument, including ``Body.at()`` and ``Topos.at()``.

* Users who want to specify a target directory when downloading a file
  will now create their own loader object, instead of having to specify
  a special keyword argument for every download::

    load = api.Loader('~/ephemeris-files')
    load('de421.bsp')

0.6.1
-----

* Users can now supply a target ``directory`` when downloading a file::

    load('de421.bsp', directory='~/ephemerides')

* Fix: removed inadvertent dependency on the Pandas library.

* Fix: ``load()`` was raising a ``PermissionError`` on Windows after a
  successful download when it tried to rename the new file.

0.6
---

* Skyfield now generates its own estimate for ``delta_t`` if the user
  does not supply their own ``delta_t=`` keyword when specifying a date.
  This should make altitude and azimuth angles much more precise.

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

* | Ask for planet positions with ``body.at(t)`` instead of ``body(t)``.

* Per IAU 2012 Resolution B2, Skyfield now uses lowercase *au* for the
  astronomical unit, and defines it as exactly 149 597 870 700 meters.
  While this API change is awkward for existing users, I wanted to make
  the change while Skyfield is still pre-1.0.  If this breaks a program
  that you already have running, please remember that a quick ``pip``
  ``install`` ``skyfield==0.4`` will get you up and running again until
  you have time to edit your code and turn ``AU`` into ``au``.

0.4
---

* To prevent confusion, the :meth:`~Time.astimezone()`
  and :meth:`~Time.utc_datetime()` methods
  have been changed to return only a ``datetime`` object.
  If you also need a leap second flag returned,
  call the new methods :meth:`~Time.astimezone_and_leap_second()`
  and :meth:`~Time.utc_datetime_and_leap_second()`.

0.3
---

* The floating-point values of an angle
  ``a.radians``, ``a.degrees``, and ``a.hours``
  are now attributes instead of method calls.


.. _Pandas API Reference: http://pandas.pydata.org/pandas-docs/stable/api.html
