
.. image:: _static/logo.png

.. rst-class:: motto

   *Elegant Astronomy for Python*

.. testsetup::

    from skyfield import api
    def now():
        """Return a constant "now"."""
        return api.JulianDate(utc=(2013, 9, 22, 14))
    api.now = now

.. testcode::

    from skyfield.api import earth, mars, now
    ra, dec, d = earth(now()).observe(mars).radec()
    print(ra)
    print(dec)

.. testoutput::

    9h 14m 50.35s
    17deg 13' 2.6"

.. testcode::

    boston = earth.topos('71.0636 W', '42.3583 N')
    alt, az, d = boston(now()).observe(mars).apparent().altaz()
    print(alt)
    print(az)

.. testoutput::

    64deg 44' 50.1"
    184deg 17' 16.9"

.. toctree::
   :maxdepth: 2

   introduction
   precision
