====================================
 Planets, and Choosing an Ephemeris
====================================

.. currentmodule:: skyfield.api

As you install Skyfield,
``pip`` will also download and install
the `de421 <https://pypi.python.org/pypi/de421>`_ package for Python.
It provides the DE421 Jet Propulsion Laboratory ephemeris
that lists high-precision planetary positions
for the years 1900 through 2050.
Its moderate 27 MB size, high accuracy, and use of very recent data
made it the best choice for Skyfield’s default ephemeris.

The DE421 ephemeris is used automatically
when you import the planets from the API module:

.. testcode::

    from skyfield.api import load

There are four other ephemerides (the plural of “ephemeris”)
available from the Python Package Index
that you might consider installing.
Here is the full list:

=========  ======== ========== ==============
Ephemeris    Size      Years       Issued
=========  ======== ========== ==============
DE405       52.1 MB  1600–2200 May 1997
DE406      170.0 MB -3000–3000 May 1997
DE421       13.0 MB  1900–2050 February 2008
DE422      519.6 MB -3000–3000 September 2009
DE423       34.6 MB  1800–2200 February 2010
=========  ======== ========== ==============

`DE405 <https://pypi.python.org/pypi/de405>`_
  A standard ephemeris from 1997 that has seen wide use
  and requires 54 MB to cover the years 1600–2200.

`DE406 <https://pypi.python.org/pypi/de406>`_
  A longer-term source of planetary positions, also published in 1997,
  that has enjoyed very widespread use
  and that stands behind several recent editions
  of the United States Naval Observatory’s Astronomical Almanac.
  Covers the years -3000–3000 at the cost of 190 MB of storage.

`DE421 <https://pypi.python.org/pypi/de421>`_
  The Skyfield default, discussed above,
  issued in 2008 which covers years 1900–2050.

`DE422 <https://pypi.python.org/pypi/de422>`_
  The most recent general-purpose long-period ephemeris from the JPL,
  issued in 2009 and covering the six thousand years
  from -3000–3000 but weighing in at 531 MB of storage.

`DE423 <https://pypi.python.org/pypi/de406>`_
  The most recent JPL ephemeris that has been packaged for Python,
  issued in 2010 and whose particular goal
  was high accuracy for Mercury and Venus
  since it was used for the MESSENGER mission to those planets.
  It covers the years 1800–2200.

To use one of these alternative ephemerides,
simply install it with the ``pip`` ``install`` command,
import it, and wrap it in a Skyfield :class:`~jpllib.Ephemeris`.
You can then use its planet attributes to turn dates into positions,
just like the planets that you can import
from the :mod:`~api` module:

.. testsetup::

    import numpy as np
    np.set_printoptions(suppress=True, precision=4)

.. testcode::

    from skyfield import api
    from skyfield.jpllib import Ephemeris

    de423 = api.load('de421.bsp')
    de430 = api.load('de423.bsp')
    jd = api.JulianDate(utc=(1993, 5, 15))

    print('DE421: {0}'.format(de423['mercury'].at(jd).position.km))
    print('DE423: {0}'.format(de430['mercury'].at(jd).position.km))

.. testoutput::

    DE421: [ 31387022.7906  33032548.9719  14337050.8211]
    DE423: [ 31387022.7693  33032548.7756  14337050.9898]

These positions might look identical
and make you wonder why the JPL even bothered to issue a new ephemeris.
But remember that these positions are in kilometers.
Read all the way down to the very last digits —
these predicted positions do differ by one or two hundred meters.
While a few-hundred-meter difference
might be impossible to see with your telescope,
it is quite relevant for JPL professionals
who need to know exactly where to aim delicate spacecraft for landings.

.. testcleanup::

    import numpy as np
    np.set_printoptions(precision=8)
