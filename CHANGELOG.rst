
Changelog
=========

1.26 — 2020 August 1
--------------------

* The official ∆T files on NASA’s FTP server have stopped receiving
  updates — they have no new data beyond February, the start of the
  global pandemic.  Unless they are updated by next February, older
  versions of Skyfield will unfortunately download the files all over
  again every time :meth:`~skyfield.iokit.Loader.timescale()` is called
  (unless the ``builtin=True`` parameter is provided).  To make Skyfield
  less fragile going forward:

  1. The loader’s :meth:`~skyfield.iokit.Loader.timescale()` method now
     defaults to ``builtin=True``, telling it to use the ∆T and leap
     second files that ship with Skyfield internally.  To download new
     ∆T files from NASA and the leap second file from the International
     Earth Rotation Service, specify ``builtin=False``.

  2. The concept of an “expired” file has been removed from ``load()``.
     Skyfield is now much simpler: if a file with the correct name
     exists, Skyfield uses it.  See :ref:`downloading-timescale-files`
     if you still want your application to check the age of your
     timescale files and automatically download new ones.

* The `ICRF.separation_from()` method now officially supports the
  combination of an array of positions with a single reference position!
  Its previous support for that combination was, alas, accidental, and
  was broken with the 1.23 release.
  `#414 <https://github.com/skyfielders/python-skyfield/issues/414>`_
  `#424 <https://github.com/skyfielders/python-skyfield/issues/424>`_

* A prototype :func:`~skyfield.magnitudelib.planetary_magnitude()`
  routine has been added with support for several planets.
  `#210 <https://github.com/skyfielders/python-skyfield/issues/210>`_

* The ``utc`` timezone that Skyfield returns in Python datetimes is now
  either the Python Standard Library’s own UTC object, if it supplies
  one, or else is defined by Skyfield itself.  Skyfield no longer
  silently tries importing the whole ``pytz`` package merely to use its
  UTC object — which also means that the timezone returned by Skyfield
  longer offers the non-standard ``localize()`` method.
  `#413 <https://github.com/skyfielders/python-skyfield/issues/413>`_

1.25 — 2020 July 24
-------------------

* Added :func:`~skyfield.data.stellarium.parse_constellations()`
  and :func:`~skyfield.data.stellarium.parse_star_names()`
  to load Stellarium star names and constellation lines.
  Constellation lines are featured in a new example script
  :ref:`neowise-chart` that produces a finder chart
  for comet C/2020 F3 NEOWISE.

* The Hipparcos star catalog should now load faster, having switched
  behind the scenes to a higher performance Pandas import routine.

* Fixed the ability of :meth:`~skyfield.timelib.Timescale.utc()` to
  accept a Python ``datetime.date`` object as its argument.
  `#409 <https://github.com/skyfielders/python-skyfield/issues/409>`_

* Slightly lowered the precision of two tests when they detect that
  Python is compiled for a 32-bit processor, so the test suite can
  succeed when contributors package Skyfield for 32-bit Linux.
  `#411 <https://github.com/skyfielders/python-skyfield/issues/411>`_

1.24 — 2020 July 20
-------------------

* Added methods :meth:`~skyfield.timelib.Timescale.from_datetime()` and
  :meth:`~skyfield.timelib.Timescale.from_datetimes()` to the
  :class:`~skyfield.timelib.Timescale` class, to better advertise the
  ability to build a Skyfield time from a Python ``datetime`` — an ability
  that was previously overloaded into the ``year`` parameter of the
  :meth:`~skyfield.timelib.Timescale.utc()` method (where it is still
  supported for backwards compatibility, but no longer documented).

* Fix: improved the accuracy with which velocity is converted between
  the Earth-fixed ITRF frame that rotates with the Earth and the
  inertial GCRS frame that does not.  In particular, this should make
  Earth satellite velocities more accurate.

1.23 — 2020 July 9
------------------

* Added :doc:`kepler-orbits` support
  for generating the positions of comets and asteroids
  from Minor Planet Center data files.

* Added :func:`~skyfield.positionlib.ICRF.is_behind_earth()` to
  determine whether a celestial object is blocked from an Earth
  satellite’s view by the Earth itself.

* Replaced the awkward and hard-to-explain ``rough_period`` search
  parameter with the conceptually simpler ``step_days`` parameter, and
  updated the instructions in :doc:`searches` to match.

* Made the :meth:`~skyfield.iokit.Loader.tle_file()` import method less
  strict about Earth satellite names: any text on the line before two
  lines of TLE data is now saved as the satellite name.  A parameter
  ``skip_names=True`` turns this off if, for particular TLE files, this
  leads to unwanted text being saved.

1.22 — 2020 Jun 8
-----------------

* Skyfield’s improved time precision (stored internally as two floats)
  is now used in computing ephemeris positions, Earth orientation, and
  light-travel time, producing position angles which change much more
  smoothly over time on a sub-milliarcsecond scale.

* :doc:`searches` is now documented for custom events that users define
  themselves, instead of only being documented for the official
  pre-written :doc:`almanac` functions.  Not only discrete events but
  also maxima and minima are now officially supported and documented,
  thanks to a rewrite of the underlying code.

* Time objects no longer cache the nutation and precession matrices,
  since they are never used again after being multiplied together to
  create the equinox-of-date rotation matrix.  This should save 144
  bytes for each time in a :class:`~skyfield.timelib.Time` array.

* It is now possible to :ref:`from-satrec` thanks to a new Earth
  satellite constructor method.
  `#384 <https://github.com/skyfielders/python-skyfield/issues/384>`_

* Added :meth:`~skyfield.iokit.Loader.build_url()` that returns the URL
  from which Skyfield will download a file.
  `#382 <https://github.com/skyfielders/python-skyfield/issues/382>`_

* Added :meth:`~skyfield.jpllib.SpiceKernel.close()` to support
  applications that need to do fine-grained resource management or whose
  testing framework check for dangling open files.
  `#374 <https://github.com/skyfielders/python-skyfield/issues/374>`_

* Skyfield’s dependency list now asks for “jplephem” version 2.13 or
  later.  Skyfield 1.21, alas, could incur a ``Module not found`` error
  when importing ``jplephem.exceptions`` if a user had an old “jplephem”
  version already installed.
  `#386 <https://github.com/skyfielders/python-skyfield/issues/386>`_

1.21 — 2020 May 29
------------------

* Added :func:`~skyfield.positionlib.ICRF.is_sunlit()` to determine
  whether Earth satellites in orbit are in Earth’s shadow or not, thanks
  to a pull request from Jesse Coffey.

* Added :func:`~skyfield.positionlib.position_of_radec()`
  to replace the poorly designed ``position_from_radec()``.

* Skyfield :class:`~skyfield.timelib.Time` objects now have microsecond
  internal accuracy, so round trips to and from Python datetimes should
  now preserve all the microsecond digits.

* The :meth:`~skyfield.timelib.Time.utc_strftime()` method now rounds to
  the nearest minute or second if it sees that either minutes or seconds
  are the smallest unit of time in the format string.

* The 6 numbers in the sequence ``t.utc`` can now be accessed by the
  attribute names ``year``, ``month``, ``day``, ``hour``, ``minute``,
  and ``second``.

* Nutation routines should now be faster and have a smaller memory
  footprint, thanks to a rewrite that uses more optimized NumPy calls.
  `#373 <https://github.com/skyfielders/python-skyfield/issues/373>`_

* Thanks to Jérôme Deuchnord, the exception raised when asking for a
  position out-of-range of a JPL ephemeris now shows the calendar dates
  for which the ephemeris is valid and carries several useful attributes.
  `#356 <https://github.com/skyfielders/python-skyfield/pull/356>`_

1.20 — 2020 April 24
--------------------

* Erik Tollerud contributed a fix for a deprecation warning about SSL
  from the most recent versions of Python (“cafile, cpath and cadefault
  are deprecated, use a custom context instead”).  The file download
  routine now auto-detects which mechanism your Python supports.
  `#363 <https://github.com/skyfielders/python-skyfield/pull/363>`_

* Added an ``elevation_m`` argument to
  :meth:`~skyfield.planetarylib.PlanetaryConstants.build_latlon_degrees()`.

1.19 — 2020 April 23
--------------------

* To hopefully fix the ``SSL: CERTIFICATE_VERIFY_FAILED`` errors that
  some users encounter when downloading timescale files, Skyfield has
  taken the risk of switching away from your system’s SSL certificates
  to the certificate bundle from the ``certifi`` package.
  `#317 <https://github.com/skyfielders/python-skyfield/issues/317>`_

* Added a new almanac routine for finding :ref:`lunar-nodes`.
  `#361 <https://github.com/skyfielders/python-skyfield/issues/361>`_

* Gave topos objects a new :meth:`~skyfield.toposlib.Topos.itrf_xyz()`
  method that returns their raw ITRF coordinates.
  `#354 <https://github.com/skyfielders/python-skyfield/issues/354>`_

* Fixed the sign of the velocity vector when two vectors are directly
  geometrically subtracted.
  `#355 <https://github.com/skyfielders/python-skyfield/issues/355>`_

1.18 — 2020 March 26
--------------------

* Deprecated the old hybrid-key satellite dictionary returned by
  ``load.tle()`` in favor of a simple list returned by the new
  :meth:`~skyfield.iokit.Loader.tle_file()` routine.
  `#345 <https://github.com/skyfielders/python-skyfield/issues/345>`_

* The almanac :func:`~skyfield.searchlib.find_discrete()` routine no
  longer returns extraneous values in its second return value if no
  changes of state were found.
  `#339 <https://github.com/skyfielders/python-skyfield/issues/339>`_
  `#351 <https://github.com/skyfielders/python-skyfield/issues/351>`_

* Added documentation and support for computing lunar libration.
  `#80 <https://github.com/skyfielders/python-skyfield/issues/80>`_

1.17 — 2020 February 2
----------------------

* Upgraded to a new version of the ``sgp4`` Python library that, when
  possible, uses the fast official C++ implementation of SGP4.

* Added a :meth:`~skyfield.sgp4lib.EarthSatellite.find_events()` Earth
  satellite method that finds the times at which a satellite rises,
  culminates, and sets.

* Improved the logic behind the :doc:`almanac` routines to avoid rare
  situations in which a cluster of nearly identical times would be
  produced for what should really be considered a single event.
  `#333 <https://github.com/skyfielders/python-skyfield/issues/333>`_

* Fixed the :meth:`~skyfield.timelib.Time.utc_strftime()` method so it
  does not report that every day in all of recorded history is a Monday.
  `#335 <https://github.com/skyfielders/python-skyfield/issues/335>`_

1.16 — 2019 December 20
-----------------------

* Added basic :doc:`planetary` support, enough to compute the position
  of a given latitude and longitude on the surface of the Moon.
  `#79 <https://github.com/skyfielders/python-skyfield/issues/79>`_
  `#124 <https://github.com/skyfielders/python-skyfield/issues/124>`_
  `#258 <https://github.com/skyfielders/python-skyfield/issues/258>`_

* Added :func:`~skyfield.almanac.oppositions_conjunctions()` for finding
  the dates when a planet is at opposition and conjunction with the sun.

* Added :func:`~skyfield.trigonometry.position_angle_of()` for computing
  astronomical position angles.

1.15 — 2019 November 20
-----------------------

* Changed the URL for the Hipparcos catalog, because the VizieR archives
  FTP server is no longer responding.
  `#301 <https://github.com/skyfielders/python-skyfield/issues/301>`_

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

* Added a :func:`~skyfield.almanac.risings_and_settings()` function for
  computing rising and setting times.
  `#271 <https://github.com/skyfielders/python-skyfield/issues/271>`_

1.13 — 2019 October 10
----------------------

* Provided a constellation lookup routine through
  :func:`~skyfield.api.load_constellation_map()`.

* Added a ``position_from_radec()`` function.

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

* Fix: :func:`~skyfield.searchlib.find_discrete()` was generating empty
  arrays of search dates, upsetting the astronomy code, if the start and
  end dates were very close together.
  `#240 <https://github.com/skyfielders/python-skyfield/issues/240>`_

1.10 — 2019 February 2
----------------------

* Fix: teach Skyfield the new format of the Naval Observatory ∆T data
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

* Both of the loader methods :meth:`~skyfield.iokit.Loader.open()` and
  ``tle()`` now accept not just URLs but also plain local file paths;
  they correctly re-download a remote file if “reload=True” is
  specified; and they allow specifying a different local “filename=”
  than the one at the end of the URL.

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
  load method downloads ∆T and leap second data from official data
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

* To prevent confusion, the :meth:`~skyfield.timelib.Time.astimezone()`
  and :meth:`~skyfield.timelib.Time.utc_datetime()` methods
  have been changed to return only a ``datetime`` object.
  If you also need a leap second flag returned,
  call the new methods
  :meth:`~skyfield.timelib.Time.astimezone_and_leap_second()`
  and :meth:`~skyfield.timelib.Time.utc_datetime_and_leap_second()`.

0.3
---

* The floating-point values of an angle
  ``a.radians``, ``a.degrees``, and ``a.hours``
  are now attributes instead of method calls.


.. _Pandas API Reference: http://pandas.pydata.org/pandas-docs/stable/api.html
