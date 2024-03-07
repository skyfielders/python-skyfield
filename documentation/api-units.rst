=======================
 API Reference — Units
=======================

When you ask positions to return distances, velocities, or angles,
they return instances of the following classes.

.. testsetup::

    from skyfield.units import Angle, Distance

.. currentmodule:: skyfield.units

.. autoclass:: Distance
   :members:

.. autoclass:: Velocity
   :members:

.. autoclass:: Angle
   :members:

.. autoclass:: AngleRate
   :members:

.. autoclass:: Rate
   :members:

.. _Formatting angles:

-----------------
Formatting angles
-----------------

To display an angle as decimal degrees or hours,
ask the angle for its ``.hours`` or ``.degrees`` attribute
and then use any normal Python mechanism for formatting a float.
For example:

.. testcode::

    ra, dec = Angle(hours=5.5877286), Angle(degrees=-5.38731536)

    print('RA {:.8f} hours'.format(ra.hours))
    print('Dec {:+.8f} degrees'.format(dec.degrees))

.. testoutput::

    RA 5.58772860 hours
    Dec -5.38731536 degrees

If you let Skyfield do the formatting instead,
then hours are split into 60 minutes of 60 seconds each,
and degrees are split into 60 arcminutes of 60 arcseconds each:

.. testcode::

    print('RA', ra)
    print('Dec', dec)

.. testoutput::

    RA 05h 35m 15.82s
    Dec -05deg 23' 14.3"

If you want more control over the display of minutes and seconds,
you can call an angle’s “hours as a string” method :meth:`~Angle.hstr`
or “degrees as a string” method :meth:`~Angle.dstr`.
The simplest adjustment you can make
is to specify the number of decimal ``places``
that will be shown in the seconds field.

.. testcode::

    print('RA', ra.hstr(places=4))
    print('Dec', dec.dstr(places=4))

.. testoutput::

    RA 05h 35m 15.8230s
    Dec -05deg 23' 14.3353"

In each of these examples
you can see that Skyfield marks arcminutes
with the ASCII apostrophe ``'``
and arcseconds
with the ASCII quotation mark ``"``.
Using plain ASCII lets Skyfield
support as many operating systems and output media as possible.
But it would be more correct
to denote arcseconds and arcminutes
with the Unicode symbols PRIME and DOUBLE PRIME,
and to use the Unicode DEGREE SIGN to mark the whole number:

−5°23′14.3″

If you want to override Skyfield’s default notation
to create either the string above, or any other notation,
then give :meth:`~Angle.hstr` or :meth:`~Angle.dstr`
a ``format=`` string of your own.
It should use the syntax of Python’s
`str.format() <https://docs.python.org/3/library/string.html#formatstrings>`_
method.
For example,
here’s the exact string you would use
to format an angle in degrees, arcminutes, and arcseconds
using the traditional typographic symbols discussed above:

.. testcode::

    print(dec.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″'))

.. testoutput::

    -5°23′14.3″

(Note that the leading ``u``, for “Unicode”, is only mandatory
in Python 2, not Python 3.)

Skyfield will call your string’s format method
with these six arguments:

{0}
    An ASCII hyphen ``'-'`` if the angle is negative,
    else the empty string.
    If you want positive angles to be decorated with a plus sign,
    try using ``{0:+>1}`` which tells Python,
    “display positional parameter 0,
    padding the field to one character wide
    if it’s less than one character already,
    and use the ``+`` character to do the padding.”

{1}
    The number of whole hours or degrees.

{2}
    The number of whole minutes.

{3}
    The number of whole seconds.

{4}
    The fractions of a second.
    Be sure to pad this field
    to the number of ``places=`` you’ve requested,
    or else a fraction like ``.0012`` will format incorrectly as ``.12``.
    If you have asked for ``places=3``, for example,
    you’ll want to display this field as ``{4:03}``.
    (See also the next item.)

{5}
    The number of ``places=`` you requested,
    which you will probably use like ``{4:0{5}}``
    when formatting field 4.
    You can use this
    in case you might not know the number of places ahead of time,
    for example if the number of places is configured by your user.

It would have been nice if the ``format=`` string
were the first option to :meth:`~Angle.hstr` or :meth:`~Angle.dstr`
so that its keyword name could be omitted
but, alas, it was only added in Skyfield 1.39,
by which point the other options had already grabbed the first spots.
