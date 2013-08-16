
===================
 Development Notes
===================

These short topics introduce new developers to some of the core concepts
of the Skyfield library, to help them be productive and to make good use
of the conventions laid down in the existing code base.

Scalar Theory
-------------

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

1. A normal Python `float` or `int` is a scalar, and lacks all of the
   special attributes that NumPy arrays enjoy.

   >>> n = 3.4
   >>> type(n)
   <type 'float'>
   >>> hasattr(n, 'shape')
   False

2. But NumPy also offers its own form of scalar: objects that have the
   same attributes and methods as arrays, but that have the empty tuple
   `()` as their `.shape` instead of specifying one or more dimensions.

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
