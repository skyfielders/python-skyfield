
===============================
 API Reference — Opening Files
===============================

.. currentmodule:: skyfield.iokit

See :doc:`files` for an explanation of how Skyfield programs
use an instance of the `Loader` class described below
to download and open the data files they need in order to operate.

The loader class
================

.. autoclass:: Loader
   :members:

   .. attribute:: directory

      The directory where this loader looks when trying to open a file,
      and where it downloads files that have not been downloaded yet.

Standalone functions
====================

.. autofunction:: load_file

.. autofunction:: parse_tle_file
