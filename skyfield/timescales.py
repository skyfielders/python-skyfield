from numpy import array, einsum, rollaxis, sin
from .constants import T0, S_DAY
from .framelib import ICRS_to_J2000 as B
from .nutationlib import compute_nutation
from .precessionlib import compute_precession

# Much of the following code is adapted from the USNO's "novas.c".

_sequence = ['tdb', 'tt', 'ut1', 'utc']
_sequence_indexes = {name: i for i, name in enumerate(_sequence)}

extra_documentation = """

        This function takes a Julian date as argument. Either provide a
        `jd=` keyword argument supplying a `JulianDate` if you already
        have one, or supply a float or array of floats with one of the
        following keyword arguments:

        `tdb=` - Barycentric Dynamical Time
        `tt=`  - Terrestrial Time
        `ut1=` - Universal Time
        `utc=` - Coordinated Universal Time

"""

def takes_julian_date(function):
    """Wrap `function` so it accepts the standard Julian date arguments.

    A function that takes two arguments, `self` and `jd`, may be wrapped
    with this decorator if it wants to support optional auto-creation of
    its `jd` argument by accepting all of the same keyword arguments
    that the JulianDate constructor itself supports.

    """
    def wrapper(self, jd=None, tdb=None, tt=None, ut1=None, utc=None,
                delta_t=0.0):
        if not isinstance(jd, JulianDate):
            jd = JulianDate(jd=jd, tdb=tdb, tt=tt, ut1=ut1, utc=utc,
                            delta_t=delta_t)
        else:
            pass  # TODO: verify that they provided a JulianDate instance
        return function(self, jd)
    wrapper.__name__ = function.__name__
    synopsis, blank_line, description = function.__doc__.partition('\n\n')
    wrapper.__doc__ = ''.join(synopsis + extra_documentation + description)
    return wrapper

def _wrap(a):
    if hasattr(a, 'shape'):
        return a
    if hasattr(a, '__len__'):
        return array(a)
    return array((a,))

def _convert(a):
    if not hasattr(a, 'shape') and hasattr(a, '__len__'):
        a = array(a)
    return a, getattr(a, 'shape', ())

class JulianDate(object):
    """Julian date.

    Attributes:

    `tdb` - Barycentric Dynamical Time
    `tt`  - Terrestrial Time
    `ut1` - Universal Time
    `utc` - Coordinated Universal Time

    """
    def __init__(self, jd=None, tdb=None, tt=None, ut1=None,
                 utc=None, delta_t=0.0):

        self.delta_t, ignored_shape = _convert(delta_t)

        if jd is not None:
            # TODO: we should actually presume that datetime's are utc
            #assert isinstance(jd, datetime)
            dt = jd
            tdb = julian_date(dt.year, dt.month, dt.day, dt.hour,
                              dt.minute, dt.second + dt.microsecond * 1e-6)
        if tdb is not None:
            self.tdb, self.shape = _convert(tdb)
        if tt is not None:
            self.tt, self.shape = _convert(tt)
        if ut1 is not None:
            self.ut1, self.shape = _convert(ut1)
        if utc is not None:
            self.utc, self.shape = _convert(utc)

        if not self.__dict__:
            raise ValueError('you must supply either tdb= tt= ut1= or utc=')

    def __getattr__(self, name):

        # Cache of several expensive functions of time.

        if name == 'P':
            self.P = P = compute_precession(self.tdb)
            return P

        if name == 'PT':
            self.PT = PT = rollaxis(self.P, 1)
            return PT

        if name == 'N':
            self.N = N = compute_nutation(self)
            return N

        if name == 'NT':
            self.NT = NT = rollaxis(self.N, 1)
            return NT

        if name == 'M':
            self.M = M = einsum('ij...,jk...,kl...->il...', self.NT, self.PT, B.T)
            return M

        if name == 'MT':
            self.MT = MT = rollaxis(self.M, 1)
            return MT

        # Conversion between timescales.

        delta_t = self.delta_t
        d = self.__dict__
        i = _sequence_indexes.get(name, None)
        if i is None:
            raise AttributeError('no such attribute {!r}'.format(name))

        _TDB, _TT, _UT1, _UTC = (0, 1, 2, 3)

        tdb = d.get('tdb')
        tt = d.get('tt')
        ut1 = d.get('ut1')
        utc = d.get('utc')  # wrong! TT -(const)-> TAI -(leap secs)-> UTC?

        if i >= _TT:
            if (tt is None) and (tdb is not None):
                self.tt = tt = tdb - tdb_minus_tt(tdb) * S_DAY
                if i == _TT:
                    return tt

            if i >= _UT1:
                if (ut1 is None) and (tt is not None):
                    self.ut1 = ut1 = tt - delta_t / S_DAY
                    if i == _UT1:
                        return ut1

                if i == _UTC:
                    if (ut1 is not None):
                        die
                        utc = self.utc = 6.7
                        return utc

        if (ut1 is None) and (utc is not None):
            die
            ut1 = self.ut1 = 3.4
            if i == _UT1:
                return ut1

        if (tt is None) and (ut1 is not None):
            self.tt = tt = ut1 + self.delta_t * S_DAY
            if i == _TT:
                return tt

        self.tdb = tdb = tt + tdb_minus_tt(tt) * S_DAY
        return tdb

def julian_date(year, month=1, day=1, hour=0.0, minute=0.0, second=0.0):
    janfeb = month < 3
    return ((second / 60.0 + minute) / 60.0 + hour) / 24.0 - 0.5 + (
            day - 32075
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            )

def cal_date(jd):
    """Convert Julian Day `jd` into a Gregorian year, month, day, and hour."""
    jd = jd + 0.5

    hour = jd % 1.0 * 24.0
    k = int(jd) + 68569
    n = 4 * k // 146097;

    k = k - (146097 * n + 3) // 4
    m = 4000 * (k + 1) // 1461001
    k = k - 1461 * m // 4 + 31
    month = 80 * k // 2447
    day = k - 2447 * month // 80
    k = month // 11

    month = month + 2 - 12 * k
    year = 100 * (n - 49) + m + k

    return year, month, day, hour

def tdb_minus_tt(jd_tdb):
    """Computes TT corresponding to a TDB Julian date."""

    t = (jd_tdb - T0) / 36525.0

    # Expression given in USNO Circular 179, eq. 2.6.

    return (0.001657 * sin ( 628.3076 * t + 6.2401)
          + 0.000022 * sin ( 575.3385 * t + 4.2970)
          + 0.000014 * sin (1256.6152 * t + 6.1969)
          + 0.000005 * sin ( 606.9777 * t + 4.0212)
          + 0.000005 * sin (  52.9691 * t + 0.4444)
          + 0.000002 * sin (  21.3299 * t + 5.5431)
          + 0.000010 * t * sin ( 628.3076 * t + 4.2490))
