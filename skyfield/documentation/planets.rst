
=============
 The Planets
=============

As you install Skyfield,
``pip`` will also download and install
the `de421 <https://pypi.python.org/pypi/de421>`_ Python package.
It provides the DE421 Jet Propulsion Laboratory ephemeris
that lists high-precision planetary positions
for the years 1900 through 2050.
Its moderate 27 MB size and use of very recent data
seemed to make it the best choice for Skyfield’s default ephemeris.

It is DE421 predictions that you will receive
when you simply import the planets from the API module:

.. testcode::

    # These use the JPL DE421 ephemeris

    from skyfield.api import sun, earth, moon  # etc

There are four other ephemerides (the plural of “ephemeris”)
available from the Python Package Index
that you might consider installing.

`DE405 <https://pypi.python.org/pypi/de405>`_
  A standard ephemeris from 1997 that has seen wide use
  and covers the years 1600 though 2200.

`DE406 <https://pypi.python.org/pypi/de406>`_
  A longer-term source of planetary positions, also published in 1997,
  that has enjoyed very widespread use
  and that stands behind several recent editions
  of the United States Naval Observatory’s Astronomical Almanac.
  Covers the years -3000 through 3000 at the cost of 190 MB of storage.

`DE421 <https://pypi.python.org/pypi/de421>`_
  The Skyfield default, discussed above,
  covering years 1900 through 2050.

`DE422 <https://pypi.python.org/pypi/de422>`_
  The most recent general-purpose long-period ephemeris from the JPL,
  covering the six thousand years
  from -3000 to 3000 but weighing in at 531 MB of storage.

`DE423 <https://pypi.python.org/pypi/de406>`_
  The most recent JPL ephemeris that has been packaged for Python,
  whose particular goal was high accuracy for Mercury and Venus
  since it was used for the MESSENGER mission to those planets.

To use one of these alternative ephemerides
simply install it with the ``pip`` ``install`` command,
import it, and wrap it in a Skyfield :class:`~jpllib.Ephemeris`.
You can then use its planet attributes to turn dates into positions,
just like the planets that you can import
from the :mod:`~skyfield.api` module:

.. testcode::

    import de423
    from skyfield.jpllib import Ephemeris

    eph = Ephemeris(de423)
    print eph.earth(utc=(1993, 5, 15)).position.AU

.. testoutput::

    
