
==============
 Design Notes
==============

.. currentmodule:: skyfield.units

This document collects various notes
about the API design of the Skyfield library,
as well as the code samples that illustrate them
and make sure that they keep working.

Position classes with coordinate methods
========================================

One of the biggest sources of confusion in PyEphem was that it offered
only a single routine to generate right ascension and declination
coordinates, whose arguments selected both what kind of position to
compute — astrometric or apparent — and also what coordinate system to
use to represent those different positions in the sky.

Users were chronically confused about which position they were asking
for, and what kind of coordinates were being used to name that position.

Skyfield therefore keeps these concepts strictly separate, by following
these rules in its API design:

* Positions are a big deal.  After all, in real life, two different
  positions are two different places up in the sky.  If beneath the
  night sky you point with your arm at one position, and then at another
  position that is not equal to the first, you will have to physically
  move your arm from the first to the second.

* Coordinates are not as big a deal.  They are merely names.  Given a
  position, you can name it with ICRS coordinates or equinox-of-date
  coordinates; with equatorial or ecliptic or galactic coordinates; but
  even as you step through all of those choices of coordinate, rattling
  off different numbers for your audience, your arm will remain in the
  exact same position.  Your finger will be pointing at exactly one
  place in the sky that all those different coordinates designate.

Therefore Skyfield deems each position as substantial enough to deserve
its own Python object.  If you convert an astrometric position into an
apparent position, you will be returned a new object to represent that
different place in the sky::

    apparent = astrometric.apparent()

But coordinates are mere names and so do not earn a separate object for
each choice of coordinate.  Instead, each kind of coordinate is simply a
different method call on the position you have generated::

    astrometric.radec()
    astrometric.ecliptic_latlon()
    astrometric.galactic_latlon()

This approach stands in contrast to both PyEphem and its one method that
generated all combinations of position and coordinate, and with AstroPy,
which dignifies different coordinates with different objects even if
those two coordinates designate exactly the same place in the sky.
Skyfield instead tries to use the object-versus-method distinction in
Python to provide extra signal to the amateur about when a program is
talking about a truly different location in the sky, versus when it is
merely generating different numbers to designate a single location.

Angle Classes
=============

Years of experience with maintaining PyEphem suggested that many users
who are first dabbling in astronomy do not realize that right ascension
is not typically expressed in degrees, and that users will frequently
become confused — generating a long support tail — if the API makes it
too natural to express right ascension using degrees.

Therefore, when faced with a request for degrees or hours, the
:class:`Angle` class tries to protect the user from accessing the wrong
one unless they have indicated that it is the value they really want.
By default, all angles assume that degrees are the proper unit, and have
no problem with the user accessing them as degrees:

>>> from skyfield.units import Angle
>>> a = Angle(degrees=-54.0)
>>> a.degrees
-54.0
>>> print(a.dstr())
-54deg 00' 00.0"
>>> a.dms()
(-54.0, -0.0, -0.0)
>>> a.signed_dms()
(-1.0, 54.0, 0.0, 0.0)

But exceptions will trigger if we try to express the angle in hours:

>>> a.hours
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in degrees, not hours; if you want to use hours anyway, then please use the attribute _hours
>>> a.hstr()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in degrees, not hours; if you want to use hours anyway, then call signed_hms() with warn=False
>>> a.hms()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in degrees, not hours; if you want to use hours anyway, then call signed_hms() with warn=False
>>> a.signed_hms()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in degrees, not hours; if you want to use hours anyway, then call signed_hms() with warn=False

The remedies suggested by these exceptions indeed work:

>>> print(a._hours)
-3.6
>>> print(a.hstr(warn=False))
-03h 36m 00.00s
>>> a.hms(warn=False)
(-3.0, -36.0, -0.0)
>>> a.signed_hms(warn=False)
(-1.0, 3.0, 36.0, 0.0)

An angle has the opposite behavior if it is given a preference of hours
when it is created.  In that case, fetching its value in hours values
becomes the easy operation:

>>> a = Angle(degrees=-54.0, preference='hours')
>>> print(a.hours)
-3.6
>>> print(a.hstr())
-03h 36m 00.00s
>>> a.hms()
(-3.0, -36.0, -0.0)
>>> a.signed_hms()
(-1.0, 3.0, 36.0, 0.0)

And it is now degrees that take a bit of extra effort to retrieve.

>>> a.degrees
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in hours, not degrees; if you want to use degrees anyway, then please use the attribute _degrees
>>> a.dstr()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in hours, not degrees; if you want to use degrees anyway, then call dstr() with warn=False
>>> a.dms()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in hours, not degrees; if you want to use hours anyway, then call signed_hms() with warn=False
>>> a.signed_dms()
Traceback (most recent call last):
  ...
WrongUnitError: this angle is usually expressed in hours, not degrees; if you want to use hours anyway, then call signed_hms() with warn=False

>>> a._degrees
-54.0
>>> print(a.dstr(warn=False))
-54deg 00' 00.0"
>>> a.dms(warn=False)
(-54.0, -0.0, -0.0)
>>> a.signed_dms(warn=False)
(-1.0, 54.0, 0.0, 0.0)

Importing NumPy
===============

Code examples in the Skyfield documentation that need NumPy
will always import it as ``np``
as that is the standard practice in the wider SciPy community.
Examples will then use NumPy features
through fully qualified names like ``np.array``
since that is how users —
especially users new to the scientific Python ecosystem —
should be advised to structure their own code.

However, because Skyfield code itself
is presumed to always use NumPy
in preference to built-in Python numerics,
the hundreds of ``np.`` prefixes would add only noise.
As a consequence, Skyfield’s modules themselves simply do a
``from`` ``numpy`` ``import`` of any names that they need.

Skyfield strives to support old versions of both Python and of NumPy,
because many users in industry and government
cannot upgrade their supporting libraries whenever they want.
So the unit tests in CI are run against NumPy 1.11.3
to maintain compatibility with that older version.

Rotation matrices or state transformation matrices?
===================================================

Instead of keeping position and velocity in separate 3-vectors of
floats, the SPICE library from the JPL concatenates them both into a
single 6-vector.  It can then express the transform of one reference
frame into another by a single 6×6 matrix.  This is clever, because a
transformed velocity is the sum of both the frame’s native velocity and
also a contribution from the angle through which the position vector is
swept.  This is very cleverly captured by the 6×6 matrix; the comments
in ``frmchg`` illustrate its structure::

           -               -
          |                 |
          |    R        0   |
          |                 |
          |                 |
          |   dR            |
          |   --        R   |
          |   dt            |
          |                 |
           -               -

The top rows say “the position is simply rotated, with no contribution
from the velocity,” while the bottom rows say “the velocity is rotated,
then added to the position × the rate the frame is rotating.”

Since an aggregate frame transform can then be constructed by simply
multiplying a series of these 6×6 matrices, a temptation arises: if
Skyfield frame objects adopted the same convention, they would only have
to carry a single transformation matrix.

The answer is no.  Skyfield does not use this technique.

To understand why, observe the waste that happens when using the above
matrix: fully one-quarter of the multiplies and something like one-half
the adds always create zeros.  The SPICE system corrects this by using a
one-off implementation of matrix multiplication ``zzmsxf`` that in fact
treats the operation as smaller 3×3 operations.  Its comments note::

        -            -    -            -
       |      |       |  |      |       |
       |   R2 |   0   |  |   R1 |   0   |
       |      |       |  |      |       |
       | -----+------ |  | -----+------ |  =
       |      |       |  |      |       |
       |   D2 |   R2  |  |   D1 |   R1  |
       |      |       |  |      |       |
        -            -    -            -

        -                              -
       |                  |             |
       |   R2*R1          |     0       |
       |                  |             |
       | -----------------+------------ |
       |                  |             |
       |   D2*R1 + R2*D1  |   R2*R1     |
       |                  |             |
        -                              -

If the cost of efficiency is the additional cost and complication of
breaking down the 6×6 so as to discard one-quarter of it and do pairwise
operations between the remaining three quarters, then Skyfield chooses
to not perform the aggregation in the first place.

So in Skyfield let’s keep the matrices ``R`` and ``dR/dt`` in the first
diagram always separate.  Then we can perform the exact 3×3 operations
that SPICE does but without what in Skyfield would be a disaggregation
step beforehand plus an aggregation step after.

Cross products
==============

How can we compute a cross product while remaining agnostic about
whether the two vectors we have been handed have a second dimension?

>>> from numpy import array, cross
>>> a = array([1, 2, 3])
>>> b = array([6, 4, 5])

>>> p = array([7, 8, 9])
>>> q = array([1, 0, 2])

The simple cross product between our 3-vectors:

>>> cross(a, p, 0, 0).T
array([-6, 12, -6])
>>> cross(b, q, 0, 0).T
array([ 8, -7, -4])

Now we combine our vectors into stacks of values across the final
dimension of a matrix:

>>> ab = array([a, b]).T
>>> ab
array([[1, 6],
       [2, 4],
       [3, 5]])
>>> pq = array([p, q]).T
>>> pq
array([[7, 1],
       [8, 0],
       [9, 2]])
>>> cross(ab, pq, 0, 0).T
array([[-6,  8],
       [12, -7],
       [-6, -4]])

(Wondering how to do fundamental_arguments() without reshape:)

>>> from numpy import array, cross, matrix, float64
>>> m = array([4.0,5.0,6.0])
>>> arg1 = float64(2.0)
>>> argn = array([2.0, 3.0])
>>> arg1 * m.reshape(3,)
array([ 8., 10., 12.])
>>> argn * m.reshape(3, 1)
array([[ 8., 12.],
       [10., 15.],
       [12., 18.]])
