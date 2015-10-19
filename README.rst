
====================================
 Welcome to the Skyfield Repository
====================================

Skyfield is a pure-Python astronomy package
that is compatible with both Python 2 and 3
and makes it easy to generate high precision research-grade
positions for planets and Earth satellites.

::

   from skyfield.api import load, now

   planets = load('de421.bsp')
   earth, mars = planets['earth'], planets['mars']

   jd = now()
   position = earth.at(jd).observe(mars)
   ra, dec, distance = position.radec()

   print(ra)
   print(dec)
   print(distance)

::

   10h 47m 56.24s
   +09deg 03' 23.1"
   2.33251 au

Its only binary dependency is NumPy.
Skyfield can usually be installed with::

    pip install skyfield

Here are the essential project links:

* `Home page and documentation
  <http://rhodesmill.org/skyfield>`_.

* `Contributing to Skyfield
  <https://github.com/skyfielders/python-skyfield/blob/master/Contrib.rst>`_.

* `Skyfield package <https://pypi.python.org/pypi/skyfield>`_
  on the Python Package Index.

* `Issue tracker
  <https://github.com/brandon-rhodes/python-skyfield/issues>`_
  on GitHub.
