
===============
 Kepler Orbits
===============

Skyfield now offers basic support for computing the position
of a comet or minor planet
whose elliptical, parabolic, or hyperbolic orbit
is provided as Kepler orbital elements.

Beware that the internal routines supporting Kepler orbits
are rudimentary and subject to change —
only the interface documented here
is guaranteed to work in future Skyfield versions.

Skyfield loads orbital elements from text files using the Pandas library.
Install it before trying any of the the examples below::

    pip install pandas

Comets
======

The `IAU Minor Planet Center <https://www.minorplanetcenter.net/>`_
distributes a ``CometEls.txt`` file
of orbital elements for predicting comet positions.
The file is plain text,
so feel free to open it with a text editor
to see the comets for which it offers orbital elements.
To build a dataframe of comets:

.. testsetup::

    def len(x):
        return 864

.. testcode::

    from skyfield.api import load
    from skyfield.data import mpc

    with load.open(mpc.COMET_URL) as f:
        comets = mpc.load_comets_dataframe(f)

    print(len(comets), 'comets loaded')

.. testoutput::

    864 comets loaded

Since the comets file has no explicit expiration date,
``load.open()`` will only download the file once.
Subsequent calls re-open the copy of the file already on your filesystem.
To force a fresh download and receive updated orbits and new comets,
pass ``reload=True``.

To generate a comet’s position,
first select its row from dataframe.
There are several Pandas techniques for selecting rows,
but most Skyfield users will simply index their dataframe by comet designation.

.. testcode::

    # Index by designation for fast lookup.
    comets = comets.set_index('designation', drop=False)

    # Sample lookups.
    row = comets.loc['1P/Halley']
    row = comets.loc['C/1995 O1 (Hale-Bopp)']

When computing the position of a comet from Earth,
there is a complication:
cometary orbits are not measured from the Solar System barycenter
but are instead centered on the Sun.
You will therefore need to add the *barycenter→Sun* vector
to the *Sun→comet* vector
to produce a position that you can pass to the ``observe()`` method,
which always measures positions from the Solar System barycenter.

.. testcode::

    # Generating a position.

    from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN

    ts = load.timescale()
    eph = load('de421.bsp')
    sun, earth = eph['sun'], eph['earth']

    comet = sun + mpc.comet_orbit(row, ts, GM_SUN)

    t = ts.utc(2020, 5, 31)
    ra, dec, distance = earth.at(t).observe(comet).radec()
    print(ra)
    print(dec)

.. testoutput::

    23h 59m 16.85s
    -84deg 46' 57.8"

Hopefully Skyfield will in the future support generating positions
for whole arrays of comets in a single efficient operation,
but for now your code should expect to operate on one comet at a time.

Minor Planets
=============

There are nearly a million minor planets
in the `IAU Minor Planet Center <https://www.minorplanetcenter.net/>`_’s
database of orbital elements,
thanks to the prodigious output of automated sky surveys
over the past few decades.

The database can be downloaded as a single ``MPCORB`` —
“Minor Planet Center orbits” —
file that offers each minor planet’s orbital elements as plain text.
But the raw file requires a bit of preprocessing
before Skyfield is ready to load it:

* The first 43 lines of the file are paragraphs that explain its contents,
  state the terms under which software programs may include the data,
  and provide links to documentation.
  Skyfield will need these non-data lines ignored or removed.

* While an uncompressed version of the file is available for download,
  most users opt to download the 55 MB compressed version
  from the Minor Planet Center
  to save bandwidth and storage.
  Decompressing the full 190 MB of data stored inside
  can require more than 1 second of computing time
  depending on your platform and processing speed.

* The complete catalog lists nearly 1 million objects,
  which can take several seconds to load and index.

For all of these reasons,
it usually makes the most sense to download, uncompress, and filter the file
before starting your application.

If your operating system provides tools for pattern matching,
they might be the fastest tool for selecting specific orbits.
Here’s how to extract the orbits
for the first four asteroids to be discovered —
(1) Ceres, (2) Pallas, (3) Juno, and (4) Vesta —
on a Linux system::

    zgrep -P '^(00001|00002|00003|00004) ' MPCORB.DAT.gz > MPCORB.excerpt.DAT

If your operating system lacks such tools,
you can build them yourself using Python.
Note that mass operations that Python implements in C,
like reading an entire file’s contents with ``read()``
and scanning the full contents with a regular expression ``findall()``,
will be much faster than using a Python loop to read every line.
Here’s an example script for performing the same search
as the ``zgrep`` command shown above:

.. include:: ../../design/mpc_make_excerpt.py
   :literal:

The same four asteroid orbits could then be extracted with::

    python mpc_make_excerpt.py 00001 00002 00003 00004 > MPCORB.excerpt.DAT

Note that the minor planets file has no explicit expiration date,
so ``load.open()`` in the above script
will only download the file once.
Subsequent calls re-open the copy of the file already on your filesystem.
To force a fresh download, pass ``reload=True``.

In either case, the resulting file — shorn of its text header,
and containing only minor planet orbits —
is ready for Skyfield to load.

.. testcode::

    with load.open('MPCORB.excerpt.DAT') as f:
        minor_planets = mpc.load_mpcorb_dataframe(f)

    print(minor_planets.shape[0], 'minor planets loaded')

.. testoutput::

    4 minor planets loaded

As was also demonstrated in the previous section on comets,
you can ask Pandas to index the dataframe by minor planet designation
for quick lookup.

.. testcode::

    # Index by designation for fast lookup.
    minor_planets = minor_planets.set_index('designation', drop=False)

    # Sample lookups.
    row = minor_planets.loc['(1) Ceres']
    row = minor_planets.loc['(4) Vesta']

Finally,
generating a position involves the same maneuver necessary for comets:
since minor planet orbits are centered on the Sun,
the Sun’s position vector must be added to theirs
to build a complete position.

.. testcode::

    ceres = sun + mpc.mpcorb_orbit(row, ts, GM_SUN)
    ra, dec, distance = earth.at(t).observe(ceres).radec()
    print(ra)
    print(dec)

.. testoutput::

    05h 51m 45.85s
    +22deg 38' 50.2"
