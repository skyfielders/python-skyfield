
==============
 Design Notes
==============

.. currentmodule:: skyfield.units

Angle Classes
-------------

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
>>> a.dstr()
'-54deg 00\' 00.0"'
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

>>> a._hours
-3.6
>>> a.hstr(warn=False)
'-03h 36m 00.00s'
>>> a.hms(warn=False)
(-3.0, -36.0, -0.0)
>>> a.signed_hms(warn=False)
(-1.0, 3.0, 36.0, 0.0)

An angle has the opposite behavior if it is given a preference of hours
when it is created.  In that case, fetching its value in hours values
becomes the easy operation:

>>> a = Angle(degrees=-54.0, preference='hours')
>>> a.hours
-3.6
>>> a.hstr()
'-03h 36m 00.00s'
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
>>> a.dstr(warn=False)
'-54deg 00\' 00.0"'
>>> a.dms(warn=False)
(-54.0, -0.0, -0.0)
>>> a.signed_dms(warn=False)
(-1.0, 54.0, 0.0, 0.0)
