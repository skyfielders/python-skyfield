
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

Comets
======

The Minor Planet Center distributes a ``CometEls.txt`` file
of orbital elements for predicting comet positions.
The file is in plain text,
so feel free to open it with a text editor
to see the comets for which it offers orbital elements.
Skyfield can import it if you first install the Pandas library::

    pip install pandas

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

To generate a comet’s position,
first select its row from dataframe.
There are several Pandas techniques for selecting rows,
but most Skyfield users will simply index their dataframe by comet name.

.. testcode::

    # Index by name for fast lookup.
    comets = comets.set_index('name', drop=False)

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

    from skyfield.keplerlib import KeplerOrbit

    ts = load.timescale(builtin=True)
    eph = load('de421.bsp')
    comet = eph['sun'] + KeplerOrbit.from_comet_row(ts, row)

    t = ts.utc(2020, 5, 31)
    ra, dec, distance = eph['earth'].at(t).observe(comet).radec()
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


