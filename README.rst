
====================================
 Welcome to the Skyfield Repository
====================================

Skyfield is a pure-Python astronomy package
that is compatible with both Python 2 and 3
and makes it easy to generate high precision research-grade
positions for planets and Earth satellites.

.. code-block:: python

   from skyfield.api import load

   planets = load('de421.bsp')
   earth, mars = planets['earth'], planets['mars']

   ts = load.timescale()
   t = ts.now()
   position = earth.at(t).observe(mars)
   ra, dec, distance = position.radec()

   print(ra)
   print(dec)
   print(distance)

The result:

::

   10h 47m 56.24s
   +09deg 03' 23.1"
   2.33251 au

Skyfield’s only binary dependency is NumPy.
Once that is available, Skyfield can usually be installed with::

    pip install skyfield

Here are the essential project links:

* `Home page and documentation
  <http://rhodesmill.org/skyfield>`_.

* `Installing Skyfield
  <http://rhodesmill.org/skyfield/installation.html>`_.

* `Contributing to Skyfield
  <https://github.com/skyfielders/python-skyfield/blob/master/Contrib.rst>`_.

* `Skyfield package <https://pypi.python.org/pypi/skyfield>`_
  on the Python Package Index.

* `Issue tracker
  <https://github.com/brandon-rhodes/python-skyfield/issues>`_
  on GitHub.
