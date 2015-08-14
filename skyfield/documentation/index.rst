
==========
 Skyfield
==========

.. rst-class:: motto

   *Elegant Astronomy for Python*

.. testsetup::

    import os
    os.chdir('../..')  # same directory as de430.bsp, hopefully
    from skyfield import api
    def now():
        """Return a constant "now"."""
        return api.JulianDate(utc=(2013, 9, 22, 14))
    api.now = now

.. testcode::

    from skyfield.api import load, now
    e = load('de430.bsp')
    ra, dec, distance = e('earth').observe('mars barycenter').at(now()).radec()

    print(ra)
    print(dec)
    print(distance)

.. testoutput::

    09h 14m 50.35s
    +17deg 13' 02.6"
    2.18573 au

Skyfield is a pure-Python astronomy package
that makes it easy to generate high precision research-grade
positions for planets and Earth satellites.
You can compute either geocentric coordinates
as shown in the example above,
or topocentric coordinates specific to your location
on the Earth’s surface:

.. testcode::

    boston = e('earth').topos('42.3583 N', '71.0636 W')
    alt, az, d = boston.observe('mars barycenter').at(now()).apparent().altaz()

    print(alt)
    print(az)

.. testoutput::

    64deg 44' 50.1"
    184deg 17' 16.9"

The official documentation is available through the links
in the Table of Contents below.
You can also visit:

* Official `Python Package Index <https://pypi.python.org/pypi/skyfield>`_
  entry

* GitHub `project page <https://github.com/brandon-rhodes/python-skyfield/>`_

* GitHub `issue tracker <https://github.com/brandon-rhodes/python-skyfield/issues>`_

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   installation
   positions
   time
   planets
   stars
   earth-satellites
   api
   design
