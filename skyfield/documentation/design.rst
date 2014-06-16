
==============
 Design Notes
==============

.. currentmodule:: skyfield.units

Angle Classes
-------------

The :class:`Angle` class could have been very simple, and have provided
symmetrical hours and degrees methods that were equally easy to call.
But years of experience with maintaining PyEphem suggested that many
users who are first dabbling in astronomy do not realize that right
ascension is not typically expressed in degrees, and that users will
frequently become confused — generating a long support tail — if the API
makes it too natural to express right ascension using degrees.

So both the hours methods on angles (like :meth:`~BaseAngle.hms()` and
:meth:`~BaseAngle.hstr()`) as well as the degrees methods (like
:meth:`~BaseAngle.dms()` and :meth:`~BaseAngle.dstr()` are methods of
:class:`BaseAngle` whose two subclasses, :class:`Angle` and
:class:`HourAngle`, are instead very opinionated about how they are
expressed.  :class:`Angle` hides the hour methods so that users are
guided towards using degrees, while :class:`HourAngle` hides the degree
methods instead so that it can safely be used for right ascension.

Both subclasses do allow the opposite kind of measure to be expressed,
but the user has to say something like ``dms_anyway()`` instead of
simply asking for ``dms()`` if the angle would not normally be expressed
as degrees.  It is possible that the string ``anyway`` should be
replaced with some more clear indication that this is an override that
violates typical astronomical convention.  An early version of the API
had ``dammit`` instead — which, it was thought, would properly express
the frustration of the user who wanted degrees and who had just been
told that they could not get them with a simple ``dms()`` call — but it
was feared that the word would not make its users’ Python code look
professional.

One last issue was that the ``repr()`` of right ascension angle objects
looked misleading:

>>> from skyfield.units import HourAngle
>>> HourAngle(hours=12.5)
HourAngle(hours=(12.0, 30.0, 0.0))

This looks wrong because the class name :class:`HourAngle`, which was
simply chosen to indicate “an angle that prefers to be printed out using
hours instead of degrees,” sounds like the technical term “hour angle”
which is a different physical measurement than a right ascension.

It is possible that :class:`HourAngle` should simply be renamed.  Maybe
someone will suggest another name in the future.  But for now, the
solution is a feature-identical subclass whose name is what you would
actually expect for printing a right ascension.

>>> from skyfield.units import RightAscension
>>> RightAscension(hours=12.5)
RightAscension(hours=(12.0, 30.0, 0.0))
