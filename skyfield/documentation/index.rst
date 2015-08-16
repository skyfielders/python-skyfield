
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
        return api.JulianDate(utc=(2015, 8, 15, 18))
    api.now = now

.. testcode::

    from skyfield.api import load, now

    planets = load('de421.bsp')
    position = planets('earth').observe('venus').at(now())
    ra, dec, distance = position.radec()

    print(ra)
    print(dec)
    print(distance)

.. testoutput::

    09h 29m 05.03s
    +06deg 36' 01.7"
    0.288446 au

Skyfield is a pure-Python astronomy package
that makes it easy to generate high precision research-grade
positions for planets and Earth satellites.
You can compute either geocentric coordinates
as shown in the example above,
or topocentric coordinates specific to your location
on the Earth’s surface:

.. testcode::

    boston = planets('earth').topos('42.3583 N', '71.0636 W')
    position = boston.observe('venus').at(now())
    alt, az, d = position.apparent().altaz()

    print(alt)
    print(az)

.. testoutput::

    49deg 48' 40.0"
    212deg 45' 09.5"

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
