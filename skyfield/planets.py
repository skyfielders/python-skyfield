import jplephem
from numpy import array

from skyfield import earthlib, timescales
from skyfield.coordinates import ICRS

DAY_S = 24.0 * 60.0 * 60.0
KM_AU = 1.0 / earthlib.AU_KM

T0 = timescales.T0

class Planet(object):
    def __init__(self, ephemeris, jplephemeris, jplname):
        self.ephemeris = ephemeris
        self.jplephemeris = jplephemeris
        self.jplname = jplname

    def __repr__(self):
        return '<Planet %s>' % (self.jplname,)

    def __call__(self, jd):
        e = self.jplephemeris
        c = e.compute
        if self.jplname == 'earth':
            pv = c('earthmoon', jd) - c('moon', jd) * e.earth_share
        elif self.jplname == 'moon':
            pv = c('earthmoon', jd) + c('moon', jd) * e.moon_share
        else:
            pv = c(self.jplname, jd)
        pv *= KM_AU
        i = ICRS(pv[:3], pv[3:], jd)
        i.ephemeris = self.ephemeris
        return i

    def observe_from(self, observer):
        return observer.observe(self)

class Ephemeris(object):

    def __init__(self, module=None):

        if module is None:
            try:
                import de421
            except ImportError:
                raise ValueError(
                    'if you want to instantiate Ephemeris() without '
                    'providing an argument, then you must install the '
                    'default ephemeris DE421 with the command: '
                    '"pip install de421"')
            else:
                module = de421

        self.jplephemeris = jplephem.Ephemeris(module)

        self.sun = Planet(self, self.jplephemeris, 'sun')
        self.mercury = Planet(self, self.jplephemeris, 'mercury')
        self.venus = Planet(self, self.jplephemeris, 'venus')
        self.earth = Planet(self, self.jplephemeris, 'earth')
        self.moon = Planet(self, self.jplephemeris, 'moon')
        self.mars = Planet(self, self.jplephemeris, 'mars')
        self.jupiter = Planet(self, self.jplephemeris, 'jupiter')
        self.saturn = Planet(self, self.jplephemeris, 'saturn')
        self.uranus = Planet(self, self.jplephemeris, 'uranus')
        self.neptune = Planet(self, self.jplephemeris, 'neptune')
        self.pluto = Planet(self, self.jplephemeris, 'pluto')

    def compute(self, name, jd):
        return getattr(self, name)(jd)
