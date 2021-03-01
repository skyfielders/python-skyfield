
==============================================
 Planets and their moons: JPL ephemeris files
==============================================

.. Big list of files is at ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp

For planets and their moons,
NASA’s Jet Propulsion Laboratory (JPL)
offers high accuracy tables of positions
for time spans ranging from decades to centuries.
Each table is called an *ephemeris*,
from the ancient Greek word for *daily*.
Here’s how to download and open the JPL ephemeris DE421
and ask for the position of Mars:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup()

.. testcode::

    from skyfield.api import load

    ts = load.timescale()
    t = ts.utc(2021, 2, 26, 15, 19)

    planets = load('de421.bsp')  # ephemeris DE421
    mars = planets['Mars Barycenter']
    barycentric = mars.at(t)

Or you can compute the position of Mars
from another vantage point in the Solar System:

.. testcode::

    earth = planets['Earth']
    astrometric = earth.at(t).observe(mars)

For details:

* See :doc:`files` to learn more about the ``load()`` routine
  and how you can choose the directory
  to which it downloads files.
  Note that ephemeris files never receive updates,
  so once you have a file like ``de421.bsp`` on disk
  you never need to download it again.

* See :doc:`positions` to learn what you can do
  with the ``barycentric`` and ``astrometric`` positions
  computed in the code above.

.. TODO “go see PLACE to learn more about vector functions”?

To learn more about ephemeris files themselves, keep reading here.

Choosing an Ephemeris
=====================

Here are the most popular general-purpose ephemeris files,
from the JPL’s famous “DE” Development Ephemeris series.

==== =============================== ====== ================================ ======
Year Short-term ephemeris            Size   Long-term ephemeris              Size
==== =============================== ====== ================================ ======
1997 ``de405.bsp`` (1600 to 2200)    63 MB  ``de406.bsp`` (−3000 to 3000)    287 MB
2008 ``de421.bsp`` (1900 to 2050)    17 MB  ``de422.bsp`` (−3000 to 3000)    623 MB
2013 | ``de430t.bsp`` (1550 to 2650) 128 MB ``de431t.bsp`` (–13200 to 17191) 3.5 GB
     | ``de430_1850-2150.bsp``       31 MB
2020 | ``de440.bsp`` (1550 to 2650)  114 MB ``de441.bsp`` (−13200 to 17191)  3.1 GB
     | ``de440s.bsp`` (1849 to 2150) 32 MB
==== =============================== ====== ================================ ======

.. TODO Link to a discussion of negative years.

How can you choose the right ephemeris for your project?
TODO: How many points below?

1. Does the ephemeris include the planets you need?
---------------------------------------------------

You can use ``print()`` to check whether an ephemeris
lists a specific planet or moon.
Here, for example, are the targets supported by DE421:

.. testcode::

    print(planets)

.. testoutput::

    SPICE kernel file 'de421.bsp' has 15 segments
      JD 2414864.50 - JD 2471184.50  (1899-07-28 through 2053-10-08)
          0 -> 1    SOLAR SYSTEM BARYCENTER -> MERCURY BARYCENTER
          0 -> 2    SOLAR SYSTEM BARYCENTER -> VENUS BARYCENTER
          0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
          0 -> 4    SOLAR SYSTEM BARYCENTER -> MARS BARYCENTER
          0 -> 5    SOLAR SYSTEM BARYCENTER -> JUPITER BARYCENTER
          0 -> 6    SOLAR SYSTEM BARYCENTER -> SATURN BARYCENTER
          0 -> 7    SOLAR SYSTEM BARYCENTER -> URANUS BARYCENTER
          0 -> 8    SOLAR SYSTEM BARYCENTER -> NEPTUNE BARYCENTER
          0 -> 9    SOLAR SYSTEM BARYCENTER -> PLUTO BARYCENTER
          0 -> 10   SOLAR SYSTEM BARYCENTER -> SUN
          3 -> 301  EARTH BARYCENTER -> MOON
          3 -> 399  EARTH BARYCENTER -> EARTH
          1 -> 199  MERCURY BARYCENTER -> MERCURY
          2 -> 299  VENUS BARYCENTER -> VENUS
          4 -> 499  MARS BARYCENTER -> MARS

Skyfield can generate positions
for any body that an ephemeris links to target zero,
the Solar System barycenter.
This ephemeris DE421, as you can see above,
provides a segment directly linking the Solar System barycenter
with the Sun:

.. testcode::

    sun = planets['Sun']
    print(sun)

.. testoutput::

    'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 10 SUN

By contrast,
generating a position for the Moon with DE421 requires two segments.
The first segment provides the position of the Earth-Moon center of gravity,
while the second segment provides the offset from there to the Moon.

.. testcode::

    moon = planets['Moon']
    print(moon)

.. testoutput::

    Sum of 2 vectors:
     'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
     'de421.bsp' segment 3 EARTH BARYCENTER -> 301 MOON

Note that most planets are so massive compared to their moons
that you can ignore the difference
between the planet and its system barycenter.
If you want to observe Mars or Jupiter from elsewhere in the Solar System,
just ask for the ``Mars Barycenter``
or ``Jupiter Barycenter`` position instead.
The Earth-Moon system is unusual
for featuring a satellite with so much mass —
though even in that case,
their common barycenter is always inside the Earth.
Only Pluto has a satellite so massive and so distant
that the Pluto-Charon barycenter is in space between them.

.. TODO Restore discussion of the “segments” list?

2. What dates do you need?
--------------------------

Most JPL ephemeris files are mission-specific
and usually less than two centuries in length —
DE418, for example,
was `issued in support of the 2015 New Horizons Pluto encounter
<https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de418_announcement.pdf>`_
and covers only the 150 years from 1900 to 2050.

But once or twice a decade
the JPL releases general-purpose ephemeris files
that are not tied to a specific mission.
These tend to come in pairs:
a short-term ephemeris and a long-term ephemeris
produced from the same data.
You can see several examples in the table above.
The most recent example is the short-term DE440 and the long-term DE441.

The most obvious tradeoff
between a short-term ephemeris
and its long-term counterpart is file size.

But there is also a less obvious tradeoff:
the most recent short-term ephemeris files
feature more accurate Moon positions
than the long-term ephemeris files,
because each short term ephemeris
models the effect of the Moon’s recently discovered fluid core.
The long-term ephemeris files
can’t include this detail because the correction is
“not suitable for backward integration of more than a few centuries”
(in the words of the DE430/DE431 report)
and would reduce the accuracy of positions in the far past or future.

3. How big a difference does the ephemeris make?
------------------------------------------------

You can always expect more recent JPL ephemeris files to be more accurate.
But how much more accurate are they?

It might surprise you to learn that JPL ephemeris files
don’t offer hard numbers for the outer limits of their accuracy.
Instead, the official reports —
listed below in the :ref:`Ephemeris bibliography` —
offer statements like these from the report on DE430 and DE431:

* “The orbits of the inner planets are known to subkilometer accuracy”
* “an accuracy of 0″.0002 …
  is the limiting error source for the orbits of the terrestrial planets,
  and corresponds to orbit uncertainties of a few hundred meters.”
* “The orbits of Jupiter and Saturn
  are determined to accuracies of tens of kilometers”
* “Uranus, Neptune, and Pluto … observations …
  limit position accuracies to several thousand kilometers.”

A recent paper by Nanograv Collaboration
makes the stark admission that
`“the ephemerides do not generally provide usable error representations.”
<https://arxiv.org/pdf/2001.00595.pdf>`_.

Instead, you will need to switch questions.

The key is to realize that the choice before you is never
“would I rather use the ephemeris DE430
or instead use the real position of Jupiter?”
Since you have no access to the real position,
your choice is always a choice between ephemeris files themselves:
“would I rather use DE430 or DE440?”
And this, happily, is a choice that can be quantified.
You can use Skyfield to try out both ephemeris files
and see how big the difference is.

Making an excerpt of an ephemeris
=================================

Several of the ephemeris files listed below are very large.
While most programmers will follow the example above and use DE421,
if you wish to go beyond its 150-year period
you will need a larger ephemeris.
And programmers interested in the moons of Jupiter
will need JUP310, which weighs in at nearly a gigabyte.

What if you need data from a very large ephemeris,
but don’t require its entire time span?

When you installed Skyfield another library named ``jplephem``
will have been installed.
When invoked from the command line,
it can build an excerpt of a larger ephemeris
without needing to download the entire file,
thanks to the fact that HTTP supports a ``Range:`` header
that asks for only specific bytes of a file.
For example,
let’s pull two weeks of data for Jupiter’s moons
(using a shell variable ``$u`` for the URL
only to make the command less wide here on the screen
and easier to read)::

$ u=https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup310.bsp
$ python -m jplephem excerpt 2018/1/1 2018/1/15 $u jup_excerpt.bsp

The resulting file ``jup_excerpt.bsp`` weighs in
at only 0.2 MB instead of 932 MB
but supports all of the same objects as the original JUP310
over the given two-week period::

  $ python -m jplephem spk jup_excerpt.bsp
  File type DAF/SPK and format LTL-IEEE with 13 segments:
  2458119.75..2458210.50  Jupiter Barycenter (5) -> Io (501)
  2458119.50..2458210.50  Jupiter Barycenter (5) -> Europa (502)
  2458119.00..2458210.50  Jupiter Barycenter (5) -> Ganymede (503)
  2458119.00..2458210.50  Jupiter Barycenter (5) -> Callisto (504)
  ...

You can load and use it directly off of disk
with :func:`~skyfield.iokit.load_file()`.

Closing the ephemeris file automatically
========================================

If you need to close files as you finish using them
instead of waiting until the application exits,
each Skyfield ephemeris offers a
:meth:`~skyfield.jpllib.SpiceKernel.close()` method.
You can either call it manually,
or use Python’s |closing|_ context manager
to call ``close()`` automatically when a block of code finishes:

.. |closing| replace:: ``closing()``
.. _closing: https://docs.python.org/3/library/contextlib.html#contextlib.closing

.. testcode::

    from contextlib import closing

    ts = load.timescale()
    t = ts.J2000

    with closing(planets):
        planets['venus'].at(t)  # Ephemeris can be used here

.. testcleanup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()

.. _third-party-ephemerides:

Type 1 and Type 21 ephemeris formats
====================================

If you generate an ephemeris with a tool like NASA’s
`HORIZONS <https://ssd.jpl.nasa.gov/horizons.cgi>`_ system,
it might be in a format not yet natively supported by Skyfield.
The first obstacle to opening the ephemeris
might be its lack of a recognized suffix:

.. testcode::

    load('wld23593.15')

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: Skyfield does not know how to open a file named 'wld23593.15'

A workaround for the unusual filename extension
is to open the file manually using Skyfield’s JPL ephemeris support.
The next obstacle, however, will be a lack of support
for Type 21 ephemerides in Skyfield:

.. testcode::

    from skyfield.jpllib import SpiceKernel
    kernel = SpiceKernel('wld23593.15')

.. testoutput::

    Traceback (most recent call last):
      ...
    ValueError: SPK data type 21 not yet supported

Older files with a similar format
might instead generate the complaint
“SPK data type 1 not yet supported.”

Happily, thanks to Shushi Uetsuki,
a pair of third-party libraries exist
that offer preliminary support for Type 1 and Type 21 ephemerides!

* https://pypi.org/project/spktype01/
* https://pypi.org/project/spktype21/

Their documentation already includes examples of generating raw coordinates,
but many Skyfield users will want to use them
in conjunction with standard Skyfield methods like ``observe()``.
To integrate them with the rest of Skyfield,
you will want to define a new vector function class
that calls the third-party module to generate coordinates:

.. testcode::

    from skyfield.constants import AU_KM
    from skyfield.vectorlib import VectorFunction
    from spktype21 import SPKType21

    t = ts.utc(2020, 6, 9)

    eph = load('de421.bsp')
    earth = eph['earth']

    class Type21Object(VectorFunction):
        def __init__(self, kernel, target):
            self.kernel = kernel
            self.center = 0
            self.target = target

        def _at(self, t):
            k = self.kernel
            r, v = k.compute_type21(0, self.target, t.whole, t.tdb_fraction)
            return r / AU_KM, v / AU_KM, None, None

    kernel = SPKType21.open('wld23593.15')
    chiron = Type21Object(kernel, 2002060)

    ra, dec, distance = earth.at(t).observe(chiron).radec()
    print(ra)
    print(dec)

.. testoutput::

    00h 27m 38.99s
    +05deg 57' 08.9"

Hopefully this third-party support
for Type 1 and Type 23 SPK ephemeris segments
will be sufficient for projects that need them,
until there is time for a Skyfield contributor
to integrate such support into Skyfield itself.

.. _Ephemeris bibliography:

Ephemeris bibliography
======================

.. TODO Mention python -m ... to show comment

Download directories

* For planets:

  | ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/
  | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/

* For planet moons:

  | ftp://ssd.jpl.nasa.gov/pub/eph/satellites/bsp/
  | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/

DE405 / DE406

* `JPL Planetary and Lunar Ephemerides, DE405/LE405
  <ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/de405.iom.pdf>`_
  (Standish 1998)

* `Check on JPL DE405 using modern optical observations
  <https://aas.aanda.org/articles/aas/pdf/1998/18/ds1546.pdf>`_
  (Morrison and Evans 1998)

* `CCD Positions for the Outer Planets in 1996–1997
  Determined in the Extragalactic Reference Frame
  <https://iopscience.iop.org/article/10.1086/300507/fulltext/>`_
  (Stone 1998)

* `Astrometry of Pluto and Saturn
  with the CCD meridian instruments of Bordeaux and Valinhos
  <https://www.aanda.org/articles/aa/full/2002/09/aa1965/aa1965.html>`_
  (Rapaport, Teixeira, Le Campion, Ducourant1, Camargo,
  Benevides-Soares 2002)

DE421

* `The Planetary and Lunar Ephemeris DE421
  <https://ipnpr.jpl.nasa.gov/progress_report/42-178/178C.pdf>`_
  (Folkner, Williams, Boggs 2009)

DE430 / DE431

* `The Planetary and Lunar Ephemerides DE430 and DE431
  <https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf>`_
  (Folkner, Williams, Boggs, Park, Kuchynka 2014)

* `DE430 Lunar Orbit, Physical Librations and Surface Coordinates
  <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf>`_
  (Williams, Boggs, Folkner 2013)

DE440 / DE441

* `The JPL Planetary and Lunar Ephemerides DE440 and DE441
  <https://iopscience.iop.org/article/10.3847/1538-3881/abd414>`_
  (Park, Folkner, Williams, and Boggs 2021)

Analysis mentioning several ephemerides

* `Modeling the Uncertainties of Solar-System Ephemerides
  for Robust Gravitational-Wave Searches with Pulsar Timing Arrays
  <https://arxiv.org/pdf/2001.00595.pdf>`_
  (The NANOGrav Collaboration 2020)

File format ``.bsp`` documentation

* `SPICE toolkit: SPK Required Reading
  <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>`_
  (describes ``.bsp`` files)

* `SPICE toolkit: Double Precision Array Files (DAF)
  <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>`_
  (describes binary format)
