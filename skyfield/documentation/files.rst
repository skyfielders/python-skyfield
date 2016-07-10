
==================================
 Downloading and Using Data Files
==================================

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

   This instance of the `Loader` class is initialized
   with the default options:
   it will automatically download missing files,
   

   Open a file by inspecting its file extension and loading its data
   with the correct Skyfield function or class.  If the file is not
   found in the current directory (or a custom data directory, if the
   user has build their own `Loader`) and Skyfield can guess its URL,
   then the file is downloaded first.

download and use data files that describe timescales,
the motion of planets, and the orbits of Earth satellites.
Most programs can simply use the instance of the `Loader` class
built with its default options
that is already available inside the :mod:`~skyfield.api` module:

.. testcode::

    from skyfield.api import load

But as :doc:`files` explains,
other programs will want to build their own loader
so that they have the chance to specify different options.

needs data files

look it downloads them

but it does not download them again, it reuses them
you can print the log and see

you can give it a custom data directory

you can turn off its verbosity

you can tell it to ignore expiration
