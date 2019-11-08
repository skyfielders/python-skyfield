Changelog
=========

.. currentmodule:: skyfield.positionlib

Next
----

* Added a :func:`~skyfield.almanac.dark_twilight_day()` function that
  not only handles sunrise and sunset but also all three kinds of
  twilight.
  `#225 <https://github.com/skyfielders/python-skyfield/issues/225>`_

1.14 — 2019 November 1
----------------------

* Changed the URL from which leap second files are downloaded; the
  server that previously provided them is no longer responding.
  Thanks to Richard Shaw for the pull request.
  `#296 <https://github.com/skyfielders/python-skyfield/issues/296>`_
  `#297 <https://github.com/skyfielders/python-skyfield/issues/297>`_

* Added a :func:`~skyfield.almanac.rising_setting()` function for
  computing rising and setting times.
  `#271 <https://github.com/skyfielders/python-skyfield/issues/271>`_

1.13 — 2019 October 10
----------------------

* Provided a constellation lookup routine through
  :func:`~skyfield.api.load_constellation_map()`.

* Added :func:`~skyfield.positionlib.position_from_radec()`.

* Fixed the ``apparent()`` method in the case where a single observer
  position is observing an entire vector of target positions.
  `#229 <https://github.com/skyfielders/python-skyfield/issues/229>`_

1.12 — 2019 September 2
-----------------------

* Fix: an exception was being thrown when creating a ``Loader`` pointed
  at a Windows directory for which Python’s ``os.makedirs()`` function
  returned a spurious error.
  `#283 <https://github.com/skyfielders/python-skyfield/issues/283>`_

* The internal ``reverse_terra()`` routine can now be given an
  ``iterations=0`` argument if the caller wants geocentric latitude and
  longitude.

1.11 — 2019 July 22
-------------------

* You can now call ``load.timescale(builtin=True)`` to use time scale
  files that Skyfield carries internally, instead of downloading them.
  Note that the time scale files distributed with any given version of
  Skyfield will gradually fall out of date.

* Fix: indexing a position now returns a position with an actual velocity.
  `#241 <https://github.com/skyfielders/python-skyfield/issues/241>`_

* Fix: the ``Star`` method ``from_dataframe()`` now correctly pulls
  stellar parallax data from the dataframe if available.
  `#266 <https://github.com/skyfielders/python-skyfield/issues/266>`_

* Fix: `find_discrete()` was generating empty arrays of search dates,
  upsetting the astronomy code, if the start and end dates were very
  close together.
  `#240 <https://github.com/skyfielders/python-skyfield/issues/240>`_

1.10 — 2019 February 2
----------------------

* Fix: teach Skyfield the new format of the Naval Observatory ΔT data
  file ``deltat.preds``, whose change in format caused Skyfield to start
  throwing an exception for new users.
  `#236 <https://github.com/skyfielders/python-skyfield/issues/236>`_

1.9 — 2018 September 23
-----------------------

* Added :func:`~skyfield.almanac.seasons` to the :doc:`almanac` module
  that can be used to predict solstices and equinoxes.

* Fix: the ecliptic coordinate routines no longer raise ``ValueError:
  too many values to unpack`` if they are passed a time array.
  `#207 <https://github.com/skyfielders/python-skyfield/issues/207>`_
  `#208 <https://github.com/skyfielders/python-skyfield/issues/208>`_

1.8 — 2018 September 12
-----------------------

* There is now an :doc:`almanac` module can compute the times of
  sunrise, sunset, and the phases of the moon, based on the search
  algorithms announced at my recent PyBay talk “An Import Loop and a
  Fiery Reentry.”

* Two new methods :meth:`~skyfield.positionlib.ICRF.cirs_xyz()` and
  :meth:`~skyfield.positionlib.ICRF.cirs_radec()` have been contributed
  which provide support for rotating a position into the Celestial
  Intermediate Reference System (CIRS).
  `#192 <https://github.com/skyfielders/python-skyfield/issues/192>`_

1.7 — 2018 September 3
----------------------

* Skyfield now supports loading the Hipparcos star catalog as a Pandas
  dataframe, providing the user with convenient mechanisms for looking
  up a single star by HIP number or filtering the entire catalog by
  magnitude.  See :doc:`stars` for details.

* Ecliptic coordinates can now be produced for epochs other than J2000
  thanks to a new optional parameter specifying the desired epoch for
  the :meth:`~skyfield.positionlib.ICRF.ecliptic_latlon()` method.

* A position that gives a position, velocity, and time can now be
  converted into full osculating orbital elements through the routine
  :func:`~skyfield.elementslib.osculating_elements_of()`.

* A couple of bugs in the ``load()`` routine have been fixed.
  `#193 <https://github.com/skyfielders/python-skyfield/issues/193>`_
  `#194 <https://github.com/skyfielders/python-skyfield/issues/194>`_

1.6 — 2018 July 25
------------------

* Both of the loader methods :meth:`~skyfield.iokit.Loader.load()` and
  :meth:`~skyfield.iokit.Loader.tle()` now accept not just URLs but also
  plain local file paths; they correctly re-download a remote file if
  “reload=True” is specified; and they allow specifying a different local
  “filename=” than the one at the end of the URL.

* Earth satellite objects no longer try to instantiate a timescale object
  of their own, which often kicked off an unexpected download of the three
  files needed to build a timescale.

* Satellite names are now correctly loaded from Space-Track TLE files.

* The ability to create times using Julian Dates is now better advertised,
  thanks to dedicated timescale methods whose names end in ``…_jd()``.

1.5 — 2018 July 4
-----------------

* The :meth:`~skyfield.positionlib.Geocentric.subpoint()` method
  now normalizes the longitude values it returns
  into the range −180° to 180°
  `#182 <https://github.com/skyfielders/python-skyfield/issues/182>`_
  and returns an actual elevation instead of zero.
  `#185 <https://github.com/skyfielders/python-skyfield/issues/185>`_

* Earth satellites now return a real velocity vector instead of zero.
  `#187 <https://github.com/skyfielders/python-skyfield/issues/187>`_

* Earth satellites now offer an
  :meth:`~skyfield.sgp4lib.EarthSatellite.ITRF_position_velocity_error()`
  method that returns raw ITRF coordinates for users interested in them.
  `#85 <https://github.com/skyfielders/python-skyfield/issues/85>`_

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

* You can now create a time object given the UT1 date.
  `#91 <https://github.com/skyfielders/python-skyfield/issues/91>`_

* Fractional Julian years are now available on ``Time`` objects as ``.J``.

* The parameter DUT1 is now available on ``Time`` objects as ``.dut1``.
  `#176 <https://github.com/skyfielders/python-skyfield/issues/176>`_

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
