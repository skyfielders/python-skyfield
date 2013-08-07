import jplephem
from numpy import max, min, sqrt

from skyfield import earthlib, timescales
from skyfield.coordinates import ICRS, GCRS

DAY_S = 24.0 * 60.0 * 60.0
KM_AU = 1.0 / earthlib.AU_KM
C_AUDAY = 173.1446326846693

T0 = timescales.T0

class Planet(object):
    def __init__(self, ephemeris, jplephemeris, jplname):
        self.ephemeris = ephemeris
        self.jplephemeris = jplephemeris
        self.jplname = jplname

    def __repr__(self):
        return '<Planet %s>' % (self.jplname,)

    def __call__(self, jd):
        position, velocity = self._position_and_velocity(jd.tdb)
        i = ICRS(position, velocity, jd)
        i.ephemeris = self.ephemeris
        return i

    def _position_and_velocity(self, jd_tdb):
        e = self.jplephemeris
        c = e.compute
        if self.jplname == 'earth':
            pv = c('earthmoon', jd_tdb) - c('moon', jd_tdb) * e.earth_share
        elif self.jplname == 'moon':
            pv = c('earthmoon', jd_tdb) + c('moon', jd_tdb) * e.moon_share
        else:
            pv = c(self.jplname, jd_tdb)
        pv *= KM_AU
        return pv[:3], pv[3:]

    def observe_from(self, observer):
        # TODO: should also accept another ICRS?

        jd_tdb = observer.jd.tdb
        lighttime0 = 0.0
        position, velocity = self._position_and_velocity(jd_tdb)
        vector = position - observer.position
        euclidian_distance = distance = sqrt((vector * vector).sum(axis=0))

        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            lighttime0 = lighttime
            position, velocity = self._position_and_velocity(jd_tdb - lighttime)
            vector = position - observer.position
            distance = sqrt((vector * vector).sum(axis=0))
        else:
            raise ValueError('observe_from() light-travel time'
                             ' failed to converge')

        g = GCRS(vector, velocity - observer.velocity, observer.jd)
        g.observer = observer
        g.distance = euclidian_distance
        g.lighttime = lighttime
        return g

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

    def _position_and_velocity(self, name, jd):
        return getattr(self, name)._position_and_velocity(jd)
