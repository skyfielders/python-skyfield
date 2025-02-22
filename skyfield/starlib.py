"""Python class for a distant object with, at most, proper motion."""

from numpy import array, cos, empty, isnan, outer, sin, where
from .constants import AU_KM, ASEC2RAD, C, C_AUDAY, DAY_S, T0
from .functions import length_of
from .relativity import light_time_difference
from .timelib import Time
from .units import Angle

class Star(object):
    """The position in the sky of a star or other fixed object.

    Each `Star` object specifies the position of a distant object.  You
    should provide as a right ascension and declination relative to the
    ICRS (the recent improvement upon J2000).  You can specify the
    coordinates using either floating point hours and degrees, or tuples
    that specify hour and degree fractions as minutes and seconds, or
    even full Skyfield :class:`~skyfield.units.Angle` objects (which can
    themselves be initialized using hours, degrees, or radians):

    >>> barnard = Star(ra_hours=17.963471675, dec_degrees=4.69339088889)
    >>> barnard = Star(ra_hours=(17, 57, 48.49), dec_degrees=(4, 41, 36.20))
    >>> barnard = Star(ra=Angle(hours=17.963471675),
    ...                dec=Angle(degrees=4.69339088889))

    For objects whose proper motion across the sky has been detected,
    you can supply velocities in milliarcseconds (mas) per year, and
    even a parallax and radial velocity if those are known:

    >>> barnard = Star(ra_hours=(17, 57, 48.49803),
    ...                dec_degrees=(4, 41, 36.2072),
    ...                ra_mas_per_year=-798.71,
    ...                dec_mas_per_year=+10337.77,
    ...                parallax_mas=545.4,
    ...                radial_km_per_s=-110.6)

    See `stars` for a guide to using a `Star` once you have created it.

    """
    au_km = AU_KM
    target = None

    def __init__(self, ra=None, dec=None, ra_hours=None, dec_degrees=None,
                 ra_mas_per_year=0.0, dec_mas_per_year=0.0,
                 parallax_mas=0.0, radial_km_per_s=0.0, names=(), epoch=T0):

        if ra_hours is not None:
            self.ra = Angle(hours=ra_hours)
        elif isinstance(ra, Angle):
            self.ra = ra
        else:
            raise TypeError('please provide either ra_hours=<float> or else'
                            ' ra=<skyfield.units.Angle object>')

        if dec_degrees is not None:
            self.dec = Angle(degrees=dec_degrees)
        elif isinstance(dec, Angle):
            self.dec = dec
        else:
            raise TypeError('please provide either dec_degrees=<float> or else'
                            ' dec=<skyfield.units.Angle object>')

        if isinstance(epoch, Time):
            epoch = epoch.tt
        # elif isinstance(epoch, float):
        #     pass
        # else:
        #     raise ValueError('the epoch= must be a Time object, or'
        #                      ' a floating point Barycentric Dynamical Time (TDB)')

        self.ra_mas_per_year = ra_mas_per_year
        self.dec_mas_per_year = dec_mas_per_year
        self.parallax_mas = parallax_mas
        self.radial_km_per_s = radial_km_per_s
        self.epoch = epoch
        self.names = names

        self._compute_vectors()

    def __repr__(self):
        opts = []
        for name in ['ra', 'dec',
                     'ra_mas_per_year', 'dec_mas_per_year',
                     'parallax_mas', 'radial_km_per_s',
                     'names', 'epoch']:
            value = getattr(self, name)
            if isinstance(value, Angle):
                value = value.degrees
            shape = getattr(value, 'shape', None)
            if shape:
                shapes = ','.join(str(n) for n in shape)
                opts.append('{0} shape={1}'.format(name, shapes))
            elif value:
                opts.append('{0}={1!r}'.format(name, value))
        return 'Star({0})'.format(', '.join(opts))

    @classmethod
    def from_dataframe(cls, df):
        epoch = 1721045.0 + _unwrap(df['epoch_year']) * 365.25
        return cls(
            ra_hours=_unwrap(df['ra_hours']),
            dec_degrees=_unwrap(df['dec_degrees']),
            ra_mas_per_year=_unwrap(df.get('ra_mas_per_year', 0)),
            dec_mas_per_year=_unwrap(df.get('dec_mas_per_year', 0)),
            parallax_mas=_unwrap(df.get('parallax_mas', 0)),
            epoch=epoch,
        )

    def _observe_from_bcrs(self, observer):
        position, velocity = self._position_au, self._velocity_au_per_d
        t = observer.t
        dt = light_time_difference(position, observer.xyz.au)
        if t.shape:
            position = (outer(velocity, t.tdb + dt - self.epoch).T + position).T
        else:
            position = position + velocity * (t.tdb + dt - self.epoch)
        if len(position.shape) > 1:
            if len(observer.xyz.au.shape) > 1:
                vector = position - observer.xyz.au
                vel = (observer.velocity.au_per_d.T - velocity).T
            else:
                vector = (position.T - observer.xyz.au).T
                vel = (observer.velocity.au_per_d - velocity.T).T
        else:
            vector = position - observer.xyz.au
            vel = (observer.velocity.au_per_d.T - velocity).T
        distance = length_of(vector)
        light_time = distance / C_AUDAY
        if len(position.shape) > 1 and not t.shape:
            tt = empty(position.shape[1:])
            tt.fill(t.tt)
            t = t.ts.tt_jd(tt)
        return vector, vel, t, light_time

    def _compute_vectors(self):
        """Compute the star's position as an ICRF position and velocity."""

        # Use 1 gigaparsec for stars whose parallax is zero.

        parallax = self.parallax_mas
        shape = getattr(parallax, 'shape', None)
        if shape:
            parallax = where(isnan(parallax), 0.0, parallax)
            parallax[parallax <= 0.0] = 1.0e-6
        elif parallax <= 0.0:
            parallax = 1.0e-6

        # Convert right ascension, declination, and parallax to position
        # vector in equatorial system with units of au.

        dist = 1.0 / sin(parallax * 1.0e-3 * ASEC2RAD)
        r = self.ra.radians
        d = self.dec.radians
        cra = cos(r)
        sra = sin(r)
        cdc = cos(d)
        sdc = sin(d)

        self._position_au = array((
            dist * cdc * cra,
            dist * cdc * sra,
            dist * sdc,
            ))

        # Compute Doppler factor, which accounts for change in light
        # travel time to star.

        k = 1.0 / (1.0 - self.radial_km_per_s / C * 1000.0)

        # Convert proper motion and radial velocity to orthogonal
        # components of motion with units of au/day.

        pmr = self.ra_mas_per_year / (parallax * 365.25) * k
        pmd = self.dec_mas_per_year / (parallax * 365.25) * k
        rvl = self.radial_km_per_s * DAY_S / self.au_km * k

        # Transform motion vector to equatorial system.

        self._velocity_au_per_d = array((
            - pmr * sra - pmd * sdc * cra + rvl * cdc * cra,
              pmr * cra - pmd * sdc * sra + rvl * cdc * sra,
              pmd * cdc + rvl * sdc,
              ))

def _unwrap(value):
    """Return floats untouched, but ask Series for their NumPy arrays."""
    return getattr(value, 'values', value)
