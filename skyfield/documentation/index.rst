
.. image:: _static/logo.png

.. rst-class:: motto

   *Elegant Astronomy for Python*

.. testsetup::

   import datetime as dtmodule
   class datetime(dtmodule.datetime):
       """Secret replacement for datetime."""
       @classmethod
       def now(cls):
           """Return a constant "now"."""
           return api.datetime(2013, 9, 22, 14)
   from skyfield import api
   api.datetime = datetime

.. testcode::

    from skyfield.api import datetime, earth, mars
    now = datetime.now()
    ra, dec, d = earth(utc=now).observe(mars).radec()
    print(ra)
    print(dec)

.. testoutput::

    9h 14m 50.346s
    17deg 13' 2.583"

.. testcode::

    boston = earth.topos('71.0636 W', '42.3583 N')
    alt, az, d = boston(utc=now).observe(mars).apparent().altaz()
    print(alt)
    print(az)

.. testoutput::

    64deg 44' 50.083"
    184deg 17' 16.943"

.. toctree::
   :maxdepth: 2

   introduction
   precision
