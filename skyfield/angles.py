# -*- coding: utf-8 -*-

import numpy as np
from .constants import tau

_to_degrees = 360.0 / tau
_from_degrees = tau / 360.0

_to_hours = 24.0 / tau
_from_hours = tau / 24.0

_instantiation_instructions = """to instantiate an Angle, try one of:

Angle(angle=another_angle)
Angle(radians=value)
Angle(degrees=value)
Angle(hours=value)

where `value` can be either a Python float or a NumPy array of floats"""

class BaseAngle(object):

    def __init__(self, angle=None, radians=None, degrees=None, hours=None):
        if angle is not None:
            if not isinstance(angle, BaseAngle):
                raise ValueError(_instantiation_instructions)
            self._radians = angle._radians
        elif radians is not None:
            self._radians = radians
        elif degrees is not None:
            self._radians = degrees * _from_degrees
        elif hours is not None:
            self._radians = hours * _from_hours

    def __format__(self, format_spec):
        return self.dstr()

    def radians(self):
        return self._radians

    def hours(self):
        return self._radians * _to_hours

    def hms(self):
        return sexagesimalize(self.hours())

    def hstr(self, places=2):
        sgn, h, m, s, etc = sexagesimalize(self.hours(), places)
        sign = '-' if sgn < 0.0 else ''
        return '%s%dh %dm %d.%0*ds' % (sign, h, m, s, places, etc)

    def degrees(self):
        return self._radians * _to_degrees

    def dms(self):
        return sexagesimalize(self.degrees())

    def dstr(self, places=1):
        sgn, d, m, s, etc = sexagesimalize(self.degrees(), places)
        sign = '-' if sgn < 0.0 else ''
        return '%s%ddeg %d\' %d.%0*d"' % (sign, d, m, s, places, etc)

    hours_anyway = hours
    hms_anyway = hms
    hstr_anyway = hstr

    degrees_anyway = degrees
    dms_anyway = dms
    dstr_anyway = dstr

class WrongUnitError(ValueError):

    def __init__(self, name, unit):
        usual = 'hours' if (unit == 'degrees') else 'degrees'
        self.args = ('This angle is usually expressed in {}, not {};'
                     ' if you want to express it in {} anyway, use'
                     ' {}_anyway()'.format(usual, unit, unit, name),)

class Angle(BaseAngle):

    __str__ = BaseAngle.dstr

    # Protect naive users from accidentally calling hour methods.

    def hours(self):
        raise WrongUnitError('hours', 'hours')

    def hms(self):
        raise WrongUnitError('hms', 'hours')

    def hstr(self):
        raise WrongUnitError('hstr', 'hours')

class HourAngle(BaseAngle):

    __str__ = BaseAngle.hstr

    # Protect naive users from accidentally calling degree methods.

    def degrees(self):
        raise WrongUnitError('degrees', 'degrees')

    def dms(self):
        raise WrongUnitError('dms', 'degrees')

    def dstr(self):
        raise WrongUnitError('dstr', 'degrees')

def sexagesimalize(value, places=0):
    sign = int(np.sign(value))
    value = np.absolute(value)
    power = 10 ** places
    n = int(7200 * power * value + 1) // 2
    n, fraction = divmod(n, power)
    n, seconds = divmod(n, 60)
    n, minutes = divmod(n, 60)
    return sign, n, minutes, seconds, fraction

def interpret_longitude(value):
    split = getattr(value, 'split', None)
    if split is not None:
        pieces = split()
        degrees = float(pieces[0])
        if len(pieces) > 1 and pieces[1].lower() == 'w':
            degrees = - degrees
        return degrees / 360. * tau
    else:
        return value

def interpret_latitude(value):
    split = getattr(value, 'split', None)
    if split is not None:
        pieces = split()
        degrees = float(pieces[0])
        if len(pieces) > 1 and pieces[1].lower() == 's':
            degrees = - degrees
        return degrees / 360. * tau
    else:
        return value
