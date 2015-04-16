
====================================
 Welcome to the Skyfield Repository
====================================

Skyfield is a pure-Python astronomy package
that is compatible with both Python 2 and 3
and makes it easy to generate high precision research-grade
positions for planets and Earth satellites.

::

    from skyfield.api import earth, mars, now
    ra, dec, distance = earth(now()).observe(mars).radec()

    print(ra)
    print(dec)
    print(distance)

::

    09h 14m 50.35s
    +17deg 13' 02.6"
    2.18572863461 AU

Its only binary dependency is NumPy.
Skyfield can usually be installed with::

    pip install skyfield

Here are the essential project links:

* `Home page and documentation
  <http://rhodesmill.org/skyfield>`_.

* `Contributing to Skyfield
  <https://github.com/ozialien/python-skyfield/blob/readme_collaboration/Contrib.rst>`_.

* `Skyfield package <https://pypi.python.org/pypi/skyfield>`_
  on the Python Package Index.

* `Issue tracker
  <https://github.com/brandon-rhodes/python-skyfield/issues>`_
  on GitHub.
