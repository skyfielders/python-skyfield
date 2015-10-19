====================================
 Planets, and Choosing an Ephemeris
====================================

.. currentmodule:: skyfield.api

If you are interested in observing the planets,
the Jet Propulsion Laboratory (JPL) has prepared long tables
recording where the planets were positioned in the past
and where they are predicted to be in the future.
A table of positions is called an *ephemeris*
and those supplied by the JPL are of very high accuracy.

You can download an ephemeris by giving :func:`load()` a filename.

A popular choice of ephemeris is DE421.
It is recent, has good precision,
was designed for general purpose use,
and is only 17 MB in size:

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.setup()

.. testcode::

    from skyfield.api import load
    planets = load('de421.bsp')

Once an ephemeris file has been downloaded to your current directly,
re-running your program will simply reuse the copy on disk
instead of downloading it all over again.

You can ``print()`` an ephemeris to learn which objects it supports.

.. testcode::

    print(planets)

.. testoutput::

    SPICE kernel file 'de421.bsp' has 15 segments
        JD 2414864.5 - JD 2471184.5  (1899-07-28 through 2053-10-08)
              0 -> 1    SOLAR SYSTEM BARYCENTER -> MERCURY BARYCENTER
              0 -> 2    SOLAR SYSTEM BARYCENTER -> VENUS BARYCENTER
              0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
              0 -> 4    SOLAR SYSTEM BARYCENTER -> MARS BARYCENTER
              0 -> 5    SOLAR SYSTEM BARYCENTER -> JUPITER BARYCENTER
              0 -> 6    SOLAR SYSTEM BARYCENTER -> SATURN BARYCENTER
              0 -> 7    SOLAR SYSTEM BARYCENTER -> URANUS BARYCENTER
              0 -> 8    SOLAR SYSTEM BARYCENTER -> NEPTUNE BARYCENTER
              0 -> 9    SOLAR SYSTEM BARYCENTER -> PLUTO BARYCENTER
              0 -> 10   SOLAR SYSTEM BARYCENTER -> SUN
              3 -> 301  EARTH BARYCENTER -> MOON
              3 -> 399  EARTH BARYCENTER -> EARTH
              1 -> 199  MERCURY BARYCENTER -> MERCURY
              2 -> 299  VENUS BARYCENTER -> VENUS
              4 -> 499  MARS BARYCENTER -> MARS

Bodies in JPL ephemeris files are each identified by an integer,
but Skyfield translates them so that you do not have to remember
that a code like 399 stands for the Earth and 499 for Mars.

Each ephemeris segment predicts the position of one body
with respect to another,
so several segments sometimes have to be combined
to generate a complete position.
The DE421 ephemeris shown above, for example,
can produce the position of the Sun directly.
But if you ask it for the position of Earth
then it will have to add together two distances:

* From the Solar System’s center (0) to the Earth-Moon barycenter (3)
* From the Earth-Moon barycenter (3) to the Earth itself (399)

This happens automatically behind the scenes.
All you have to say is ``planets[399]`` or ``planets['Earth']``
and Skyfield will put together a solution using the segments provided.

Here are several popular ephemerides that you can use with Skyfield.
The filenames matching ``de*``
predict the positions of many or all of the major planets,
while ``jup310.bsp`` focuses on Jupiter and its major moons:

==========  ====== ============= ==============
Ephemeris    Size      Years        Issued
==========  ====== ============= ==============
de405.bsp    63 MB  1600 to 2200 May 1997
de406.bsp   287 MB -3000 to 3000 May 1997
de421.bsp    17 MB  1900 to 2050 February 2008
de422.bsp   623 MB -3000 to 3000 September 2009
de430.bsp   128 MB  1550 to 2650 February 2010
jup310.bsp  932 MB  1900 to 2100 December 2013
==========  ====== ============= ==============

You can think of negative years, as cited in the above table,
as being almost like years BC except that they are off by one.
Historians invented our calendar back before zero was a counting number,
so AD 1 was immediately preceded by 1 BC without a year in between.
But astronomers count backwards AD 2, AD 1, 0, −1, −2, and so forth.

So if you are curious about the positions of the planets back in 44 BC,
when Julius Caesar was assassinated,
be careful to ask an astronomer about the year −43 instead.
“The fault, dear Brutus, is not in our stars, but in ourselves.”

Once you have loaded an ephemeris and used a statement like

.. testcode::

    mars = planets['Mars']

to retrieve a planet, consult the chapter :doc:`positions`
to learn about all the positions that you can use it to generate.

.. testsetup::

   __import__('skyfield.tests.fixes').tests.fixes.teardown()

