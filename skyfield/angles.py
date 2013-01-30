ASEC360 = 1296000.0
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
tau = 6.283185307179586476925287

class Angle(object):

    def __init__(self, radians):
        self.radians = radians

    def hours(self):
        return 24. / tau * self.radians

    def hms(self):
        h, ms = divmod(self.radians / tau * 24., 1.)
        m, ss = divmod(ms * 60., 1.)
        s = ss * 60.
        return h, m, s

    def degrees(self):
        return 360. / tau * self.radians

    def dms(self):
        h, ms = divmod(self.radians / tau * 360., 1.)
        m, ss = divmod(ms * 60., 1.)
        s = ss * 60.
        return h, m, s

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
