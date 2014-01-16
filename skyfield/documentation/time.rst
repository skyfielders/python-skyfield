
.. currentmodule:: skyfield.api

================
 Dates and Time
================

The foundation of Skyfield
is its specially crafted :class:`JulianDate` class,
which performs efficient conversions
among the various astronomical time scales.
It also acts as an internal cache for various date-dependent values,
so that Skyfield only has to compute each of them once.

Skyfield objects that accept date arguments give you two choices.
Once choice is to build your own :class:`JulianDate`:

.. testcode::

    from skyfield.api import JulianDate, earth

    # Building the JulianDate yourself

    jd = JulianDate(utc=(2014, 1, 1))
    print earth(jd).position.AU

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

This is always the better approach
if you are going to be using a particular date more than once,
because all of the information that the date object
caches during its first computation
will be available, for free, in all subsequent computations.

But if you are only going to use a date once,
you can skip the step of creating an object.
Instead, simply supply the arguments
that you would have used to build the date,
and Skyfield will build a :class:`JulianDate` behind the scenes
on your behalf:

.. testcode::

    # Letting Skyfield build the JulianDate

    print earth(utc=(2014, 1, 1)).position.AU

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

This results in more concise code.
The only disadvantage is that you lack a date object
to re-use in any subsequent computation.

Arrays of dates
===============

To compute a position or coordinate
over a whole range of dates,
take advantage of Skyfield’s support for NumPy arrays.
Building separate date objects is a slow process:

.. testcode::

    # Building separate dates is slow and awkward

    for day in 1, 2, 3, 4:
        jd = JulianDate(utc=(2014, 1, day))
        print earth(jd).position.AU

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]
    [-0.19179872  0.88265548  0.38254134]
    [-0.20891924  0.87936337  0.38111391]
    [-0.22597338  0.87579547  0.37956709]

If you instead construct a single :class:`JulianDate`
that holds a whole array of dates,
then explicit loops disappear from your code
and Skyfield can perform efficient vector computations
over the whole range of dates at once:

.. testcode::

    # Date arrays are concise and efficient

    days = [1, 2, 3, 4]
    jd = JulianDate(utc=(2014, 1, days))
    pos = earth(jd).position.AU
    print pos

.. testoutput::

    [[-0.17461758 -0.19179872 -0.20891924 -0.22597338]
     [ 0.88567056  0.88265548  0.87936337  0.87579547]
     [ 0.38384886  0.38254134  0.38111391  0.37956709]]

Note the shape of the resulting NumPy array.
If you unpack this array using three names,
then you get three four-element arrays
whose four values correspond
to the four dates in the :class:`JulianDate`.
These arrays are ready to be submitted to `matplotlib`_
and other scientific Python tools::

    x, y, z = pos
    plot(x, y)

If you instead slice along the second axis,
then you can retrieve an individual position for a particular date —
for example, the position corresponding to the first date
lives at index zero:

.. testcode::

    print pos[:,0]

.. testoutput::

    [-0.17461758  0.88567056  0.38384886]

If you want to check,
you can scroll up and verify that this is the same coordinate
that we generated for January 1st at the top of this document
when we were using a single date instead of an array.

.. _matplotlib: http://matplotlib.org/
