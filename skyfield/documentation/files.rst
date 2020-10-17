
==================================
 Downloading and Using Data Files
==================================

The first time you run a Skyfield program,
it will typically download one or more data files from the Internet
that provide data about planet or satellite orbits —
one file for each call the program makes to Skyfield’s ``load()`` routine.
If the program is attached to a terminal,
then a simple progress bar will be displayed
as Skyfield downloads each file.

::

   from skyfield.api import load
   planets = load('de421.bsp')
   print('Ready')

::

   [#################################] 100% de421.bsp
   Ready

The second time you run the program, however,
the program will find the data file
already sitting in the current directory.
In that case, the program will use the file on disk
without needing access to the Internet:

::

   Ready

Most programs will run just fine using the default ``load()`` routine
provided in the ``skyfield.api`` module.
But other programs may want to build their own loader
so they have the chance to override its default behaviors.

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

   from skyfield.api import Loader
   load = Loader('~/skyfield-data', verbose=False)
