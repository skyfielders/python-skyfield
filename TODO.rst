=======================
The Skyfield To-Do List
=======================

This list includes both a section on immediately-reachable goals that
are appropriate for sprints and collaboration, and longer-term goals
that the code base is not quite ready for yet but that we do not want to
forget.

Sprint Possibilities
====================

* If a user builds a ``wgs84.latlon()``, and then calls its ``at()`` to
  get a position, and then tries to use it to ``observe()``, then they
  get::

    AttributeError: 'Geocentric' object has no attribute 'observe'

  This might be what I’ve been looking for: a refutation of the idea
  that position classes should have different methods.  My idea had been
  that if only the ``Barycentric`` class had ``observe()``, that users
  couldn’t call it on the wrong sort of position.

  That plan might have worked if users all used IDEs that prevented them
  from typing ``.observe()`` on a class with no such method.  But even
  then, it would offer no help in how to get hold of a class that did
  have the method.

  But when typing in an editor, nothing stops them from typing
  ``.observe()``, and the error they get back gives zero information
  about how to fix the problem — generating a support load for the
  maintainer on GitHub.

  We should look at what it would be like for ``.observe()`` to be on
  all positions, and raise informative exceptions for positions to which
  it doesn’t apply.

* Incorporate historical polar x and y (link shared in #372):
  https://datacenter.iers.org/data/latestVersion/186_EOP_C01_2000.1846_NOW_V2013_01186.txt

* Switch to a cube-root falloff, rather than a discontinuity, for Earth
  deflection.

* Explain more about the rotation matrices used in coordinate
  transforms, as they are a source of questions on GitHub.  One user
  needed (and, happily, worked out themselves!) the dynamic reference
  frame of the ecliptic::

   ecliptic = skyfield.framelib.ecliptic_frame.rotation_at(t)

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

* Some users need to add offsets to vectors like .position and .velocity
  so should there be some official way to do so?  We should either
  document the maneuver ``.position = Distance(au=old_position.au +
  xyz)`` or (gulp!) make ``Distance`` objects susceptible of being
  modified in-place without their various units going out of sync.

* The almanac routines ``sunrise_sunset()`` and ``dark_twilight_day()``
  currently use an expensive step size — with 25 steps per day — to
  catch even brief days near the poles; see `Issue 571
  <https://github.com/skyfielders/python-skyfield/issues/571>`_ for a
  user who noticed when we tried adjusting to 6⅔ steps per day.  But
  what if instead we used the very uniform rate of the Earth’s rotation
  to zero in first on noon and midnight?  Both noon and midnight could
  surely be found within, say, one second, with a single round of
  search, right?  Then the iterative search could look for the
  transition, with zero chance of missing a, say, one-minute arctic day
  or night.

* After computing a satellite position, should the position’s `.target`
  be the satellite object itself instead of a big negative integer?
  It’s an incompatible change, but more in agreement with how other
  geocentric objects work.  Was the old behavior documented?  (Not that
  I can find!)

* One contributor says they needed to work out a routine which “Adds
  days in units of TT days, which are exactly 86400 SI seconds: (I'm
  assuming you need to pass along ts explicitly, in case you have more
  than one of them)” and would appreciate more official support::

    def offset_TT_date(ts, ti, add_days):
         return(ts.tt_jd(ti.tt + add_days))

  They also needed to “add days in units of Earth rotation.”

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
