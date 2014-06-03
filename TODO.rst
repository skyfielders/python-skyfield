=======================
The Skyfield To-Do List
=======================

This list includes both a section on immediately-reachable goals that
are appropriate for sprints and collaboration, and longer-term goals
that the code base is not quite ready for yet but that we do not want to
forget.

Sprint Possibilities
====================

* What happens if an angle that's a vector of values has .dstr() or
  .hstr() called on it?  (And, those might not be those good method
  names.)  Make it return a list of strings.

* The SciPy ``optimize`` module can find when a curve reaches zero,
  which can be used to find various interesting astronomical
  circumstances.  Examine how PyEphem uses its newton method (look in
  its ``__init__.py`) for these kinds of solutions to find out how the
  general approach works.  Write up an IPython Notebook that uses
  Skyfield together with ``optimize`` to compute as many of the
  following circumstances as possible; compare your results to the USNO
  web site to see if you are getting good answers.

  * Sunrise, local noon, and sunset.
  * New moon, full moon, and the quarters in between.
  * Solar and lunar eclipses.
  * The vernal and autumnal equinoxes.

  http://aa.usno.navy.mil/data/

* We should implement comets and asteroids using the standard formulae
  (can we find a good vector version, that will match the rest of our
  approach?) for a Keplerian orbit.

* We should add formulae for the moons of the other planets.

* There should be routines for downloading current astronomical data.
  Each routine should take an optional filename but should also have a
  default that the corresponding "load" routine uses by default.  All
  such routines, because they are likely be called later with the same
  filename as last time, should by default just skip the download if the
  file is present (and perhaps recent enough, if a ``days_old=``
  parameter is provided?) unless an ``overwrite=False`` standard
  parameter is overridden to give it the value ``True``.  (Or does the
  possibility of ``days_old=0`` give us "overwrite" for free?)  Nice
  data would be:

  * Comet orbital elements
  * Asteroid orbital elements
  * Earth-satellite orbital elements
  * (What others can we think of?)

* If anyone is a NumPy expert, I would love comments on whether my code
  is at all idiomatic, or whether I'm doing things inefficiently and in
  such a quirky way that no one else will ever understand it.

* If one of the parameters of nutationlib or precessionlib was
  incorrect, would our tests detect it?  It would be nice to have an
  automated method of determining which particular constants in the code
  base our tests are not sensitive to.


Longer-term goals
=================

* Make all objects that are `.observe()`’d from Earth include a
  sublatitude and sublongitude coordinate stating the position on Earth
  from which they appear directly overhead.  When complete, make a note
  at the PyEphem GitHub issue:

  https://github.com/brandon-rhodes/pyephem/issues/16

* When Earth Satellites are implemented, include the orbit number of a
  satellite's current position in the public attributes that are set on
  the resulting position object, as promised in PyEphem GitHub issue:

  https://github.com/brandon-rhodes/pyephem/issues/15


.. testing
     we need tests that handle both use_earth True and False.
       Similarly for other variables.
   documentation
     writing up SkyField solutions to PyEphem questions on Stack Overflow
     section on accuracy of each algorithm involved
     logo?
   performance
     Is all this vectorization worth it?
       Run a loop to compute N planet positions.
       Do the same computation using a vector of N jd's.
       Compare the runtimes under both C Python and PyPy.
       Might have to do numpypy thing; do it in skyfield/__init__.py.
       If they both show a difference, then YES it is worth it.
       Could Star() become a whole catalog of stars processed in parallel?
     What routines are taking the most time when the tests are run?
     Try to take advantage of jplephem's ability to use bundles

   Whether SGP4 passes the original library's test suite. [huh?]
