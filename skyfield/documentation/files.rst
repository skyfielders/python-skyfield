
==================================
 Downloading and Using Data Files
==================================

.. currentmodule:: skyfield.api

Your Skyfield programs will typically download two kinds of data file.

First, Skyfield will need up-to-date tables about time —
files providing recently measured values for ΔT,
future predictions of ΔT, and a table of recent leap seconds.
See :doc:`time` for an introduction to these concepts.

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

   # Depending on the order that these .rst files run in, this might
   # download the files afresh, or might find them already present:

   load = Loader('.')
   ts = load.timescale()
   planets = load('de421.bsp')

   # But this second invocation of the loader will definitely find them,
   # making the "print(load.log)" below come out right:

   load = Loader('.')
   ts = load.timescale()
   planets = load('de421.bsp')

::

   from skyfield.api import load

   ts = load.timescale()
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
provided in the :mod:`skyfield.api` module.
But other programs may want to build their own loader
so that they have the chance to specify non-default behavior.

Specifying the download directory
=================================

The default ``load()`` object saves files directly
to your current working directory —
usually the folder from which you launched your Skyfield program.
If we ask the load object that we used above
for the sequence of actions that it took,
we will see that it looked for files in the current directory:

.. testcode::

    print(load.log)

.. testoutput::

    Already exists: ./deltat.data
      Parsing with: parse_deltat_data()
      Does not expire til: 2018-01-01
    Already exists: ./deltat.preds
      Parsing with: parse_deltat_preds()
      Does not expire til: 2018-01-01
    Already exists: ./Leap_Second.dat
      Parsing with: parse_leap_seconds()
      Does not expire til: 2017-07-28
    Already exists: ./de421.bsp
      Opening with: SpiceKernel

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

.. _turning-off-downloads:

Turning off downloads for expired files
=======================================

Sometimes you want to build an application
that does not need network access in order to operate.
Half of the solution is easy: simply distribute the application
along with all of the data files that it needs,
and Skyfield will find and use the files on disk
instead of needing to download them.

But the other half of the problem is that
sometimes Skyfield will find a file on disk,
but not want to use it because the file is too old.
This can happen with each of the three time scale files,
because they become out of date after several months.

Normally, you will want Skyfield to go ahead and download new copies
so that your results are as precise as possible.
But if you think that you or your users might launch your program
when they lack access to the network,
then you can tell Skyfield to go ahead and use the files on disk
regardless of whether they are too old and have expired:

.. testcode::

   load = Loader('~/skyfield-data', expire=False)

With ``expire`` set to ``False``,
Skyfield will still try to download each file the first time
if it cannot find it in the directory the loader is using.
But on all subsequent runs, it will happily keep using those files
without ever checking whether it is time for them to expire.
