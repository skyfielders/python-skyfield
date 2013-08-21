.. Skyfield documentation master file, created by
   sphinx-quickstart on Sun Feb 17 11:09:09 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/logo.png

.. rst-class:: motto

   *Elegant Astronomy for Python*

::

   from skyfield.planets import earth, mars
   boston = earth.topos('71.0636 W', '42.3583 N')
   h = boston(jd).observe(mars).apparent().horizontal()
   print h.alt.dstr()
   print h.az.dstr()

::

   40deg 3m 49.433s
   201deg 58m 7.067s

.. toctree::
   :maxdepth: 2

   introduction
   precision
