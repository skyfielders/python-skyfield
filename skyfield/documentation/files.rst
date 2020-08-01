
==================================
 Downloading and Using Data Files
==================================

Your Skyfield programs will typically download two kinds of data file.

First, Skyfield will need up-to-date tables about time —
files providing recently measured values for ∆T,
future predictions of ∆T, and a table of recent leap seconds.
See :ref:`downloading-timescale-files` to learn more about these files.

Second, Skyfield will need data
about the objects that you want to observe:
planets, stars, or Earth satellites.
The rest of the documentation
will introduce you to each of these kinds of data source.

The result will be a program that, when run for the first time,
downloads several data files before it can perform useful operations.
If the program’s output is a terminal window,
then it will display progress bars as each file is downloaded:

.. testsetup::

   import os
   from skyfield.api import Loader

   load = Loader('.')
   ts = load.timescale(builtin=False)
   planets = load('de421.bsp')

::

   from skyfield.api import load

   ts = load.timescale(builtin=False)
   planets = load('de421.bsp')
   print('Ready')

::

   [#################################] 100% deltat.data
   [#################################] 100% deltat.preds
   [#################################] 100% Leap_Second.dat
   [#################################] 100% de421.bsp
   Ready

The second time you run the program, however,
it will find the data files already sitting in the current directory
and can start up without needing to access the network:

::

   Ready

Most programs will run just fine using the default ``load()`` object
provided in the ``skyfield.api`` module.
But other programs may want to build their own loader
so that they have the chance to specify non-default behavior.

Specifying the download directory
=================================

The default ``load()`` object saves files directly
to your current working directory —
usually the folder from which you launched your Skyfield program.

But you can instead create your own loader
that uses a different directory instead.
Simply instantiate a `Loader` with the path to the directory
where you would prefer for data files to be kept.

.. testcode::

   from skyfield.api import Loader
   load = Loader('~/skyfield-data')

Now all of your ``load()`` operations
will target that directory instead.
Note that there is no restriction
on how many `Loader` objects you can create —
feel free to create one where you put time files,
another for ephemeris files, and a third for Earth satellite TLEs,
if that makes it easier for you to keep everything organized!

Turning off the progress bar
============================

If it annoys you to have a progress bar displayed on the screen
each time that Skyfield downloads a file —
which might especially be a problem when you are using Skyfield
inside of a larger application —
you can turn the progress bars off
by building a `Loader` whose verbosity is set to false.

.. testcode::

   load = Loader('~/skyfield-data', verbose=False)
