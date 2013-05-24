# -*- coding: utf-8 -*-

ASEC360 = 1296000.0
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
tau = 6.283185307179586476925287

import numpy as np
from numpy import ndarray, array


class Angle(ndarray):

    # __array_priority__ = -1.0

    # def __init__(self, radians):
    #     self.radians = radians

    def hours(self):
        return 24. / tau * self

    def hms(self):
        return sexa(self.hours())

    def hstr(self):
        sgn, h, m, s = self.hms()
        sign = '-' if sgn < 0.0 else ''
        return '{}{}h {}m {}s'.format(sign, int(h), int(m), float(s))

    def degrees(self):
        return 360. / tau * self

    def dms(self):
        return sexa(self.degrees())

    def dstr(self):
        sgn, d, m, s = sexa(self.degrees())
        sign = '-' if sgn < 0.0 else ''
        return '{}{}deg {}m {}s'.format(sign, int(d), int(m), float(s))

    def dpretty(self):
        d, m, s = self.dms()
        return '{}°{}´{}´´'.format(int(d), int(m), float(s))


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
