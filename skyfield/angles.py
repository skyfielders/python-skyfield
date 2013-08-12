# -*- coding: utf-8 -*-
import numpy as np
from numpy import ndarray, array

from .constants import TAU

class WrongUnitError(ValueError):

    def __init__(self, name, unit):
        usual = 'hours' if (unit == 'degrees') else 'degrees'
        self.args = ('This angle is usually expressed in {}, not {};'
                     ' if you want to express it in {} anyway, use'
                     ' {}_anyway()'.format(usual, unit, unit, name),)


class BaseAngle(ndarray):

    def hours(self):
        return 24. / TAU * self

    def hms(self):
        return sexa(24. / TAU * self)

    def hours_str(self):
        sgn, w, m, s = self.hms()
        list_of_strings = self.list_repr(sgn, w, m, s, 'h')
        return list_of_strings

    def degrees(self):
        return 360. / TAU * self

    def dms(self):
        return sexa(self.degrees())

    def degrees_str(self):
        sgn, w, m, s = self.dms()
        list_of_strings = self.list_repr(sgn, w, m, s, 'deg')
        return list_of_strings

    def list_repr(self, sgn, w, m, s, type):
        sgn = np.where(sgn < 0, '-', '')
        zipped = zip(sgn, w, m, s)
        list_of_strings = []
        for angle in zipped:
            list_of_strings.append('{}{}{} {}m {:.3f}s'.format(angle[0], int(angle[1]), type, int(angle[2]), float(angle[3])))
        return list_of_strings

    def dpretty(self):
        d, m, s = self.dms()
        return '{}°{}´{}´´'.format(int(d), int(m), float(s))

    hours_anyway = hours
    hms_anyway = hms
    hstr_anyway = hours_str

    degrees_anyway = degrees
    dms_anyway = dms
    dstr_anyway = degrees_str


class Angle(BaseAngle):

    # Protect naive users from accidentally calling hour methods.

    def hours(self):
        raise WrongUnitError('hours', 'hours')

    def hms(self):
        raise WrongUnitError('hms', 'hours')

    def hours_str(self):
        raise WrongUnitError('degrees_str', 'hours')



class HourAngle(BaseAngle):

    # Protect naive users from accidentally calling degree methods.

    def degrees(self):
        raise WrongUnitError('degrees', 'degrees')

    def dms(self):
        raise WrongUnitError('dms', 'degrees')

    def degrees_str(self):
        raise WrongUnitError('degrees_str', 'degrees')


one = array([1.0])

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
        return degrees / 360. * TAU
    else:
        return value

def interpret_latitude(value):
    split = getattr(value, 'split', None)
    if split is not None:
        pieces = split()
        degrees = float(pieces[0])
        if len(pieces) > 1 and pieces[1].lower() == 's':
            degrees = - degrees
        return degrees / 360. * TAU
    else:
        return value
