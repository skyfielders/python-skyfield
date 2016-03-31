
===================================
 API Reference — Downloading files
===================================

.. currentmodule:: skyfield.iokit

Most Skyfield programs begin by importing `load`
and use it to fetch all of their data.

.. testcode::

   from skyfield.api import load

   ts = load.timescale()
   de421 = load('de421.bsp')

If the necessary files are not already in the current directory,
then the user is shown progress bars as they download. ::

   [#################################] 100% deltat.data
   [#################################] 100% deltat.preds
   [#################################] 100% Leap_Second.dat
   [#################################] 100% de421.bsp

.. function:: skyfield.iokit.load(filename)

   Open a file by inspecting its file extension and loading its data
   with the correct Skyfield function or class.  If the file is not
   found in the current directory (or a custom data directory, if the
   user has build their own `Loader`) and Skyfield can guess its URL,
   then the file is downloaded first.

.. autoclass:: Loader
   :members:
