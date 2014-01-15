
===========================
 Positions and Coordinates
===========================

Skyfield is careful to distinguish the *position* of an object
from any particular choice of *coordinates*
by which the user might want to designate that position.
There are only three basic kinds of position that Skyfield recognizes,
but several different ways in which a position
can be turned into coordinates.

Here is a quick reference to the three basic kinds of position,
together with all of the attributes and methods that they support:

.. parsed-literal::

    Three positions

    obj(time)           →  Barycentric position (BCRS)
     └─ observe(obj2)   →  Astrometric position (ICRS)
         └─ apparent()  →  Apparent position (GCRS)

    Barycentric, Astrometric, or Apparent position
     │
     ├── `position <api.html#Position.position>`_.AU         →   x, y, z
     ├── `position <api.html#Position.position>`_.km         →   x, y, z
     ├── `position.to(unit) <api.html#Distance.to>`_   →   x, y, z
     │
     ├── `velocity <api.html#Position.velocity>`_.AU_per_d   →   xdot, ydot, zdot
     ├── `velocity <api.html#Position.velocity>`_.km_per_s   →   xdot, ydot, zdot
     ├── `velocity.to(unit) <api.html#Distance.to>`_   →   xdot, ydot, zdot
     │
     ├── `radec() <api.html#Position.radec>`_             →   ra, dec, distance
     └── `radec(epoch=jd) <api.html#Position.radec>`_     →   ra, dec, distance

    Apparent position only
     │
     └── `altaz() <api.html#Position.altaz>`_             →   alt, az, distance

    Angle like ra, dec, alt, and az
     │
     ├── `radians() <api.html#Angle.radians>`_           →   6.266029488577352
     │
     ├── `hours() <api.html#Angle.hours>`_             →   23.934469599999996
     ├── `hms() <api.html#Angle.hms>`_               →   (1, 23, 56, 4, 0)
     ├── `hstr() <api.html#Angle.hstr>`_              →   '23h 56m 4.09s'
     ├── `hstr(places=4) <api.html#Angle.hstr>`_      →   '23h 56m 4.0906s'
     │
     ├── `degrees() <api.html#Angle.degrees>`_           →   359.017044
     ├── `dms() <api.html#Angle.dms>`_               →   (1, 359, 1, 1, 0)
     ├── `dstr() <api.html#Angle.dstr>`_              →   '359deg 1\' 1.4"'
     └── `dstr(places=3) <api.html#Angle.dstr>`_      →   '359deg 1\' 1.358"'

The rest of this page is simply designed to explain
all of the features outlined in the quick reference above.
All hyperlinked attributes and method names,
both in the text above and in the explanations below,
lead to the low-level :doc:`api`
which explains each option in even greater detail.

From barycentric to astrometric position
========================================

When you ask Skyfield for the position of a planet or star,
it produces a three-dimensional position
that is measured from the Solar System *barycenter* —
the center of mass around which all of the planets revolve.
It is stored as an *x*, *y*, and *z* coordinate
in the Barycentric Celestial Reference System (BCRS),
which you can view by asking Skyfield
for the :attr:`~Position.position` attribute:

.. testcode::

    from skyfield.api import earth, jupiter
    print earth(utc=(1980, 1, 1)).position.AU
    print jupiter(utc=(1980, 1, 1)).position.AU

.. testoutput::

    [-0.16287311  0.88787399  0.38473904]
    [-4.71061475  2.32932129  1.11328106]

The coordinates shown above are measured
using the Astronomical Unit (AU),
which is the average distance from the Earth to the Sun.
So the value ``-4.71`` indicates a distance
nearly five times farther from the Sun than that of the Earth.
You can, if you want, ask for these coordinates
in kilometers with the :attr:`~Position.km` attribute instead.
And if you have the third-party AstroPy package installed,
then you can convert these coordinates
into any length unit with the :meth:`~Position.to()` method.

You might think that you could determine
the position of Jupiter in the night sky
by simply subtracting these two positions
to generate the vector difference between them.
But that would ignore the fact that light takes several minutes
to travel between Jupiter and the Earth.
The image of Jupiter in our sky
does not show us where it *is*, right now,
but where it *was* — several minutes ago —
when the light now reaching our eyes or instruments
actually left its surface.

Correcting for the light-travel time
does not simply fix a minor inconvenience,
but reflects a very deep physical reality.
Not only the light from Jupiter,
but *all* of its physical effects,
arrive no faster than the speed of light!
As Jupiter tugs us with its gravity,
we do not get pulled in the direction of the “real” Jupiter —
we get tugged in the direction of its time-delayed image
hanging in the sky above us!

So Skyfield offers a :meth:`~Position.observe()` method
that carefully backdates the position of another object
to determine where it *was* when it generated the image
that we see in our sky:

.. testcode::

    astro = earth(utc=(1980, 1, 1)).observe(jupiter)
    print astro.position.AU

.. testoutput::

    [-4.54763822  1.44160883  0.72860876]

This light-delayed position is called the *astrometric* position,
and is traditionally mapped on a star chart
by the angles *right ascension* and *declination*
that you can compute using the :meth:`~Position.radec()` method
and display using their :meth:`~Angle.hstr()`
and :meth:`~Angle.dstr()` methods:

.. testcode::

    ra, dec, distance = astro.radec()
    print ra.hstr()
    print dec.dstr()
    print distance.AU

.. testoutput::

    10h 49m 38.71s
    +08deg 41' 00.6"
    4.82598384993

As we will explore in the next section,
objects never appear at exactly the position in the sky
predicted by the simple and ideal astrometric position.
But it is useful for mapping the planet
against the background of stars in a
`printed star atlas <http://www.amazon.com/s/?_encoding=UTF8&camp=1789&creative=390957&linkCode=ur2&pageMinusResults=1&suo=1389754954253&tag=letsdisthemat-20&url=search-alias%3Daps#/ref=nb_sb_noss_1?url=search-alias%3Daps&field-keywords=star%20atlas&sprefix=star+%2Caps&rh=i%3Aaps%2Ck%3Astar%20atlas&sepatfbtf=true&tc=1389754955568>`_.


