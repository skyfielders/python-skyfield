
==============
 Design Notes
==============

.. currentmodule:: skyfield.units

This document collects various notes
about the API design of the Skyfield library,
as well as the code samples that illustrate them
and make sure that they keep working.

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

Scalar Theory
=============

The Skyfield library is designed both for the use of naive users who
understand only the concept of a single date and resulting position, as
well as for users who confidently submit entire NumPy arrays of dates
for which they expect entire arrays of output.

How can a single API support both?

To begin with, all lower-level functions use pure vectors.  So when
implementing a fundamental algorithm, go ahead and assume that dates and
positions are all going to be vectors holding zero or more quantities.

It is only the top-level friendly API classes that need to support
scalars.  They each do this through two maneuvers: *wrapping* each
incoming scalar value as a length-one array, and *unwrapping* the
resulting outputs before the user sees them.

The maneuver is only slightly complicated by the fact that, when working
in a NumPy universe, there are two slightly different kinds of scalars
that the user might submit.

1. A normal Python ``float`` or ``int`` is a scalar, and lacks all of
   the special attributes that NumPy arrays enjoy.

   >>> n = 3.4
   >>> type(n).__name__
   'float'
   >>> hasattr(n, 'shape')
   False

2. But NumPy also offers its own form of scalar: objects that have the
   same attributes and methods as arrays, but that have the empty tuple
   ``()`` as their ``.shape`` instead of specifying one or more
   dimensions.

   >>> import numpy as np
   >>> m = np.float64(3.4)
   >>> m.shape
   ()

Our official approach to testing scalar-hood, therefore, is to attempt
to access the shape to determine whether it is the empty tuple, with the
empty tuple provided as the fallback if the object lacks a shape
entirely:

>>> assert getattr(n, 'shape', ()) == ()
>>> assert getattr(m, 'shape', ()) == ()
>>> assert getattr(np.array([7]), 'shape', ()) != ()

Note that, in addition, an empty array does *not* qualify as a scalar:

>>> assert getattr(np.array([]), 'shape', ()) != ()

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

