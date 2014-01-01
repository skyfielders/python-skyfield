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
        return sexa(self.hours())

    def hstr(self):
        sgn, h, m, s = sexa(self.hours())
        sign = '-' if sgn < 0.0 else ''
        return '{}{}h {}m {:.3f}s'.format(sign, int(h), int(m), float(s))

    def degrees(self):
        return self._radians * _to_degrees

    def dms(self):
        return sexa(self.degrees())

    def dstr(self):
        sgn, d, m, s = sexa(self.degrees())
        sign = '-' if sgn < 0.0 else ''
        return '{}{}deg {}m {:.3f}s'.format(sign, int(d), int(m), float(s))

    def dpretty(self):
        d, m, s = self.dms()
        return '{}°{}´{}´´'.format(int(d), int(m), float(s))

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

    # Protect naive users from accidentally calling hour methods.

    def hours(self):
        raise WrongUnitError('hours', 'hours')

    def hms(self):
        raise WrongUnitError('hms', 'hours')

    def hstr(self):
        raise WrongUnitError('hstr', 'hours')

class HourAngle(BaseAngle):

    def __format__(self, format_spec):
        return self.hstr()

    # Protect naive users from accidentally calling degree methods.

    def degrees(self):
        raise WrongUnitError('degrees', 'degrees')

    def dms(self):
        raise WrongUnitError('dms', 'degrees')

    def dstr(self):
        raise WrongUnitError('dstr', 'degrees')

def sexa(value):
    sign = np.sign(value, subok=False)
    absolute = np.absolute(value, subok=False)
    fraction, whole = np.modf(absolute)
    fraction, minutes = np.modf(fraction * 60.0)
    seconds = fraction * 60.0
    return sign, whole, minutes, seconds

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
