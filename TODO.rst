=======================
The Skyfield To-Do List
=======================

This list includes both a section on immediately-reachable goals that
are appropriate for sprints and collaboration, and longer-term goals
that the code base is not quite ready for yet but that we do not want to
forget.

Sprint Possibilities
====================

* The “finals2000A.all” file expresses DUT1 like ``0.8084178`` with 7
  digits after the decimal place.  How accurate is that?  It goes down
  to tenths of microseconds.  Roughly how much arc does the Earth rotate
  through in a tenth of a microsecond?

  (1e-7 / 24 / 60 / 60) * 360 * 60 * 60 * 1e3
  = 0.0015

  So we know the Earth’s orientation to within 0.0015 mas.  How does
  that compare to the error of storing TT in a 64-bit float?  Per the
  output of ``design/time_precision.py``, using a single float for TT
  provides steps of around 4e-5 and 5e-5 from 1950 through 2050: about
  500 times less precise than the DUT1 numbers!  Therefore, we should
  design delta-T to use both TT floats.

* Should the deflection routine really say this instead?::

            limb_angle, nadir_angle = compute_limb_angle(
                target_au + gcrs_position, gcrs_position)

* Improve the situation around “observer data”.

  * We currently special-case two components of the full
    center-to-target vector: the barycenter-to-target vector for things
    like deflection, aberration, and like magnitude being affected by
    distance from the sun; and Earth-to-target, for Earth deflection.
    First, those could reasonably live on the position object itself, as
    they are not really specific to the observer; and, second, we could
    instead simply save the components we added together.  But maybe
    that would take too much RAM?  Let’s continue to special case those
    two sub-results only, to save RAM, but from now on keep them on the
    position itself.

  * It feels like further simplifications might happen in ``VectorSum``
    but we’ll see.  Maybe ``._at()`` should finally become something
    more visible, like ``.vectors_at()``?  Its purpose is to return a
    position and velocity without building a whole position object that
    would just get thrown away in the case of a sum of two vectors, for
    which the user will only need a final position object for the
    combined sum.

* Trying to index a unit class should print help suggesting a unit be
  specified, similar to trying to iterate across one.

* Several interesting API questions arise because of
  `this Stack Overflow question <https://stackoverflow.com/questions/62654081/path-between-two-topos-locations-determine-latitude-and-longitude-where-a-giv>`_.
  Should I finally go through with renaming ``.position`` to ``.xyz``?
  Should ``.from_altaz()`` retain observer data
  so that a follow-up ``.altaz()`` returns the same coordinates?
  Should positions fully support math,
  so the computation of the Topos position can only be computed once?

* For #145, skip deflections of planets that can’t affect an observation.

* For #145, create a good syntax for combining two ephemerides.

* Expand a bit on the documentation for stars, now that they can have an
  epoch for their position.

* When I wrote `add_deflection()` and needed to know whether Jupiter
  itself is available in an ephemeris, or whether the Jupiter Barycenter
  should be used in its place, I tried writing the test `if name not in
  ephemeris:`.  But this sent the code into an infinite loop!  Why does
  `in` not work on an ephemeris object?  This should be fixed, and a
  test written to keep it fixed.

* The deflection code should really use integer identifiers instead of
  using names, which are slower because they need decoding.

* The deflection code should have a quick way to reach in and ask an
  ephemeris for a raw position in au, without having to spin up a body
  object and have it spin up a `Distance` object.

* We currently download most SPICE kernels from NAIF, but have to use
  FTP for fetching DE422.  Are the files from the two sites equivalent
  and do they have the same data?  Should we prefer one or the other?

* We should have an illustration of Earth satellite heights above the
  surface, plotted against a blue atmosphere fading out into the black
  of space as the plot goes upwards towards the top.

* The ``Timescale`` object does not currently know where its data comes
  from, so its ``repr()`` is pretty uninformative.  Should it someday be
  put in charge of loading its data from the data files specified, so
  that its ``repr()`` can print out which leap second file and delta T
  file it is using?  Should it also display how up to date the files
  are, and what leap seconds it knows about?  Also: it should be told
  the expiration date of all of its data, so that it can print it out as
  part of its ``repr()``.

* In `stars.rst`, document the other alternatives for how to set the RA
  and dec of a new Star object.

* Solar eclipses.

Adding more smarts to ephemeris handling
========================================

* When we add or subtract vectors, the new `VectorSum` needs to inherit
  an `ephemeris` from one of the segments being combined.  Right now it
  just grabs the first one.  What if, of the choices of ephemeris among
  the segments, grabbing the first one gives us an ephemeris that is
  missing several key large bodies, and so we run into an exception when
  `.apparent()` tries to compute gravitational deflection?  Should we be
  more intelligent in our choice?  Or should we even combine the various
  ephemerides our segments might offer?  Or should we specifically go
  ahead and look for the deflectors we need and try to find a segment
  for them each?

* And additionally: the error when not enough bodies are available for
  deflection maybe someday needs to be more helpful.

* Bring support for Type 1 and Type 21 JPL ephemerides into Skyfield
  from the third-party libraries where it now lives.  The routines will
  need to be vectorized and updated to use a Skyfield approach to vector
  and matrix operations.  When complete and documented, make a comment
  at: https://github.com/skyfielders/python-skyfield/issues/350

For 2.0
=======

* Remove old deprecation warnings for pre-1.0 behaviors.

* Remove support and tests for old ephemeris Python packages.
