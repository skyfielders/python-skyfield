import sys
import math
from numpy import(abs, amax, amin, arange, arccos, arctan, array, cos, cosh, 
                  cross, exp, log, nan, ndarray, newaxis, ones_like, pi, power, 
                  repeat, sin, sinh, sqrt, sum, tan, tile, zeros_like)

from skyfield.functions import length_of, dots
from skyfield.descriptorlib import reify
from skyfield.elementslib import OsculatingElements, normpi
from skyfield.units import Distance, Velocity, Angle
from skyfield.vectorlib import VectorFunction
from skyfield.data.gravitational_parameters import GM_dict
from skyfield.constants import AU_KM, DAY_S


class KeplerOrbit(VectorFunction):
    def __init__(self, 
                 position, 
                 velocity, 
                 epoch, 
                 mu_km_s=None, 
                 mu_au_d=None, 
                 center=None, 
                 target=None, 
                 center_name=None, 
                 target_name=None
        ):
        """ Calculates the position of an object using 2 body propagation
        
        Parameters
        ----------
        position : Distance
            Position vector at epoch with shape (3,)
        velocity : Velocity
            Velocity vector at epoch with shape (3,)
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_km_s : float
            Value of mu (G * M) in km^3/s^2
        mu_au_d : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for 
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        if (mu_km_s and mu_au_d) or (not mu_km_s and not mu_au_d):
            raise ValueError('Either mu_km_s or mu_au_d should be used, but not both')

        if mu_km_s:
            self._mu_km_s = mu_km_s
        elif mu_au_d:
            self._mu_km_s =  mu_au_d / AU_KM**3 * DAY_S**2
        
        self.position_at_epoch = position
        self.velocity_at_epoch = velocity
        self.epoch = epoch
        self.center = center
        self.target = target
        self.center_name = center_name
        self.target_name = target_name
        
    
    @classmethod
    def from_true_anomaly(cls, p, e, i, Om, w, v, 
                          epoch, 
                          mu_km_s=None,
                          mu_au_d=None,
                          center=None, 
                          target=None,
                          center_name=None,
                          target_name=None,
        ):
        """ Creates a `KeplerOrbit` object from elements using true anomaly
        
        Parameters
        ----------
        p : Distance
            Semi-Latus Rectum
        e : float
             Eccentricity
        i : Angle
            Inclination
        Om : Angle
            Longitude of Ascending Node
        w : Angle
            Argument of periapsis
        v : Angle
            True anomaly
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_km_s : float
            Value of mu (G * M) in km^3/s^2
        mu_au_d : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for 
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        if (mu_km_s and mu_au_d) or (not mu_km_s and not mu_au_d):
            raise ValueError('Either mu_km_s or mu_au_d should be used, but not both')
            
        if mu_au_d:
            mu_km_s = mu_au_d / AU_KM**3 * DAY_S**2
        
        position, velocity = ele_to_vec(p.km, 
                                        e, 
                                        i.radians, 
                                        Om.radians, 
                                        w.radians, 
                                        v.radians, 
                                        mu_km_s,
        )
        return cls(Distance(km=position), 
                   Velocity(km_per_s=velocity), 
                   epoch, 
                   mu_km_s, 
                   center=center, 
                   target=target,
                   center_name=center_name,
                   target_name=target_name,
        )
        

    @classmethod
    def from_mean_anomaly(cls, p, e, i, Om, w, M, 
                          epoch, 
                          mu_km_s=None,
                          mu_au_d=None,
                          center=None, 
                          target=None,
                          center_name=None,
                          target_name=None,
        ):
        """ Creates a `KeplerOrbit` object from elements using mean anomaly
        
        Parameters
        ----------
        p : Distance
            Semi-Latus Rectum
        e : float
             Eccentricity
        i : Angle
            Inclination
        Om : Angle
            Longitude of Ascending Node
        w : Angle
            Argument of periapsis
        M : Angle
            Mean anomaly
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_km_s : float
            Value of mu (G * M) in km^3/s^2
        mu_au_d : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for 
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        if (mu_km_s and mu_au_d) or (not mu_km_s and not mu_au_d):
            raise ValueError('Either mu_km_s or mu_au_d should be used, but not both')
            
        if mu_au_d:
            mu_km_s = mu_au_d  / AU_KM**3 * DAY_S**2
        
        E = eccentric_anomaly(e, M.radians)
        v = Angle(radians=true_anomaly(e, E))
        pos, vel = ele_to_vec(p.km, 
                              e, 
                              i.radians, 
                              Om.radians, 
                              w.radians, 
                              v.radians, 
                              mu_km_s,
        )
        return cls(Distance(km=pos), 
                   Velocity(km_per_s=vel), 
                   epoch, 
                   mu_km_s, 
                   center=center, 
                   target=target,
                   center_name=center_name,
                   target_name=target_name,
        )
    
    
    @classmethod
    def from_mpcorb_dataframe(cls, df, ts):
        if 'Number' in df:
            target = int(df.Number.strip('()')) + 2000000
        else:
            target = None
        
        if 'Name' not in df or df.Name == nan:
            target_name = df.Principal_desig
        else:
            target_name = df.Name
            
        p = df.a * (1 - df.e**2)
        return cls.from_mean_anomaly(p=Distance(au=p),
                                     e=df.e,
                                     i=Angle(degrees=df.i), 
                                     Om=Angle(degrees=df.Node), 
                                     w=Angle(degrees=df.Peri), 
                                     M=Angle(degrees=df.M),
                                     epoch=ts.tdb_jd(df.Epoch),
                                     mu_km_s=GM_dict[10] + GM_dict.get(target, 0),
                                     center=10,
                                     target=target,
                                     center_name='SUN',
                                     target_name=target_name,
        )
    
    
    @classmethod
    def from_comet_dataframe(cls, df, ts):
        mu_km_s = GM_dict[10]
        mu_au_d = mu_km_s / (AU_KM**3) * (DAY_S**2)
        e = df.e
        a = df.Perihelion_dist / (1 - e)
        p = a * (1 - e**2)
        n = sqrt(mu_au_d/a**3)
        peri_day = ts.tt(df.Year_of_perihelion, 0, df.Day_of_perihelion)
        epoch = ts.tt(df.Epoch_year, df.Epoch_month, df.Epoch_day)
        M = n * (epoch - peri_day)
        return cls.from_mean_anomaly(p=Distance(au=p),
                                     e=e,
                                     i=Angle(degrees=df.i),
                                     Om=Angle(degrees=df.Node), 
                                     w=Angle(degrees=df.Peri), 
                                     M=Angle(radians=M),
                                     epoch=epoch,
                                     mu_km_s=mu_km_s,
                                     center=10,
                                     # TODO: infer target SPK-ID from info in dataframe
                                     center_name='SUN',
                                     target_name=df.Designation_and_name,
        )
        
        
    def _at(self, time):
        """Propagate the KeplerOrbit to the given Time object
        
        The Time object can contain one time, or an array of times
        """
        pos, vel = propagate(self.position_at_epoch.km, 
                             self.velocity_at_epoch.km_per_s,
                             self.epoch.tt,
                             time.tt,
                             self._mu_km_s,
        )
        return pos / AU_KM, vel / AU_KM * DAY_S, None, None
        
        
    @reify
    def elements_at_epoch(self):
        return OsculatingElements(self.position_at_epoch, 
                                  self.velocity_at_epoch, 
                                  self.epoch,
                                  self._mu_km_s,
        )
    
    
    def __str__(self):
        ele = self.elements_at_epoch
        if self.target_name:
            return 'KeplerOrbit {0} {1} -> {2} {3}'.format(self.center,
                                                           self.center_name,
                                                           self.target,
                                                           self.target_name,
            )
        else:
            ele = self.elements_at_epoch
            string = 'KeplerOrbit {0} {1} -> q={2:.2}au, e={3:.1f}, i={4:.1f}°, Ω={5:.1f}°, ω={6:.1f}°'
            return string.format(self.center,
                                 self.center_name,
                                 ele.periapsis_distance.au,
                                 ele.eccentricity,
                                 ele.inclination.degrees,
                                 ele.longitude_of_ascending_node.degrees,
                                 ele.argument_of_periapsis.degrees,
            )
    
    
    def __repr__(self):
        return '<{0}>'.format(str(self))
         

def eccentric_anomaly(e, M):
    """ Iterates to solve Kepler's equation to find eccentric anomaly
    
    Based on the algorithm in section 8.10.2 of the Explanatory Supplement 
    to the Astronomical Almanac, 3rd ed. 
    """
    M = normpi(M)
    E = M + e*sin(M)
    
    max_iters = 100
    iters = 0
    while iters < max_iters:
        dM = M - (E - e*sin(E))
        dE = dM/(1 - e*cos(E))
        E = E + dE
        if abs(dE) < 1e-14: return E
        iters += 1
    else:
        raise ValueError('Failed to converge')
    

def true_anomaly(e, E):
    """Calculates true anomaly from eccentric anomaly
    
    Equation from  here step 3 here:
    https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    """
    return 2 * arctan(((1+e)/(1-e))**.5 * tan(E/2))


def ele_to_vec(p, e, i, Om, w, v, mu):
    """Calculates state vectors from orbital elements. Also checks for invalid
    sets of elements.

    Based on equations from this document:

    https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    """
    # Checks that longitude of ascending node is 0 if inclination is 0
    if isinstance(i, ndarray) or isinstance(Om, ndarray):
        if ((i==0)*(Om!=0)).any():
            raise ValueError('If inclination is 0, longitude of ascending node must be 0')
    else:
        if i==0 and Om!=0:
            raise ValueError('If inclination is 0, longitude of ascending node must be 0')

    # Checks that argument of periapsis is 0  if eccentricity is 0
    if isinstance(e, ndarray) or isinstance(w, ndarray):
        if ((e==0)*(w!=0)).any():
            raise ValueError('If eccentricity is 0, argument of periapsis must be 0')
    else:
        if e==0 and w!=0:
            raise ValueError('If eccentricity is 0, argument of periapsis must be 0')

    # Checks that true anomaly is less than arccos(-1/e) for hyperbolic orbits
    if isinstance(e, ndarray) and isinstance(v, ndarray):
        inds = (e>1)
        if (v[inds]>arccos(-1/e[inds])).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    elif isinstance(e, ndarray) and not isinstance(v, ndarray):
        inds = (e>1)
        if (v>arccos(-1/e[inds])).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    elif isinstance(v, ndarray) and not isinstance(e, ndarray):
        if e>1 and (v>arccos(-1/e)).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    else:
        if e>1 and v>arccos(-1/e):
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')

    # Checks that inclination is between 0 and pi
    if isinstance(i, ndarray):
        assert ((i>=0) * (i < pi)).all()
    else:
        assert i>=0 and i<pi
        
    r = p/(1 + e*cos(v))
    h = sqrt(p*mu)
    u = v+w

    X = r*(cos(Om)*cos(u) - sin(Om)*sin(u)*cos(i))
    Y = r*(sin(Om)*cos(u) + cos(Om)*sin(u)*cos(i))
    Z = r*(sin(i)*sin(u))

    X_dot = X*h*e/(r*p)*sin(v) - h/r*(cos(Om)*sin(u) + sin(Om)*cos(u)*cos(i))
    Y_dot = Y*h*e/(r*p)*sin(v) - h/r*(sin(Om)*sin(u) - cos(Om)*cos(u)*cos(i))
    Z_dot = Z*h*e/(r*p)*sin(v) + h/r*sin(i)*cos(u)

    # z and z_dot are independent of Om, so if Om is an array and the other
    # elements are scalars, z and z_dot need to be repeated
    if Z.size!=X.size:
        Z = repeat(Z, X.size)
        Z_dot = repeat(Z_dot, X.size)

    return array([X, Y, Z]), array([X_dot, Y_dot, Z_dot])


dpmax = sys.float_info.max

def bracket(num, end1, end2):
    num[num<end1] = end1
    num[num>end2] = end2
    return num
    
    
def find_trunc():
    denom = 2
    factr = 2
    trunc = 1
    x = 1 / denom
    while 1+x > 1:
        denom = denom * (2+factr) * (1+factr)
        factr = factr + 2
        trunc = trunc + 1
        x = 1 / denom
    return trunc


trunc = find_trunc()
odd_factorials = array([math.factorial(i) for i in range(3, trunc*2, 2)])
even_factorials = array([math.factorial(i) for i in range(2, trunc*2, 2)])


def stumpff(x):
    """Calculates Stumpff functions
    
    Based on the function toolkit/src/spicelib/stmp03.f from the SPICE toolkit, 
    which can be downloaded from naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
    """
    if (x < (-(log(2) + log(dpmax))**2)).any():
        raise ValueError('Argument below lower bound')
        
    z = sqrt(abs(x))
    
    c0 = zeros_like(x)
    c1 = zeros_like(x)
    c2 = zeros_like(x)
    c3 = zeros_like(x)
    
    low = x < -1
    c0[low] = cosh(z[low])
    c1[low] = sinh(z[low])/z[low]
    c2[low] = (1 - c0[low])/x[low]
    c3[low] = (1 - c1[low])/x[low]
    
    high = x > 1
    c0[high] = cos(z[high])
    c1[high] = sin(z[high])/z[high]
    c2[high] = (1 - c0[high])/x[high]
    c3[high] = (1 - c1[high])/x[high]
    
    mid = ~low * ~high
    n = sum(mid)
    exponents = tile(arange(0, trunc-1), [n, 1])
    odd_denominators = tile(odd_factorials, [n, 1])
    even_denominators = tile(even_factorials, [n, 1])
    numerators = repeat(x[mid][newaxis].T, trunc-1, axis=1)
    c3[mid] = sum(power(numerators, exponents)/odd_denominators, axis=1)
    c2[mid] = sum(power(numerators, exponents)/even_denominators, axis=1)
    c1[mid] = 1 - x[mid]*c3[mid]
    c0[mid] = 1 - x[mid]*c2[mid]
    
    return c0, c1, c2, c3
    

def propagate(position, velocity, t0, t1, gm):
    """Propagates a position and velocity vector with an array of times.
    
    Based on the function toolkit/src/spicelib/prop2b.f from the SPICE toolkit, 
    which can be downloaded from naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
    
    Parameters
    ----------
    position : ndarray
       Position vector with shape (3,)
    velocity : ndarray
        Velocity vector with shape (3,)
    t0 : float
        Time corresponding to `position` and `velocity`
    t1 : float or ndarray
        Time or times to propagate to
    gm : float
        Gravitational parameter in units that match the other arguments 
    """
    if gm <= 0: 
        raise ValueError("'gm' should be positive")
    if length_of(velocity) == 0:
        raise ValueError('Velocity vector has zero magnitude')
    if length_of(position) == 0:
        raise ValueError('Position vector has zero magnitude')
        
    r0 = length_of(position)
    rv = dots(position, velocity)
    
    hvec = cross(position, velocity)
    h2 = dots(hvec, hvec)
    
    if h2 == 0: 
        raise ValueError('Motion is not conical')
 
    eqvec = cross(velocity, hvec)/gm + -position/r0
    e = length_of(eqvec)
    q = h2 / (gm * (1+e))
    
    f = 1 - e
    b = sqrt(q/gm)
    
    br0 = b * r0
    b2rv = b**2 * rv
    bq = b * q
    qovr0 = q / r0
    
    
    maxc = max(abs(br0), 
               abs(b2rv),
               abs(bq),
               abs(qovr0/(bq)))
    
    if f < 0:
        fixed = log(dpmax/2) - log(maxc)
        rootf = sqrt(-f)
        logf = log(-f)
        bound = min(fixed/rootf, (fixed + 1.5*logf)/rootf)
    else:
        logbound = (log(1.5) + log(dpmax) - log(maxc)) / 3
        bound = exp(logbound)
        
    def kepler(x):
        c0, c1, c2, c3 = stumpff(f*x*x)
        return x*(br0*c1 + x*(b2rv*c2 + x*(bq*c3)))
    
    dt = t1 - t0
    
    if not isinstance(dt, ndarray):
        dt = array([dt])
        return_1d_array = True
    else:
        return_1d_array = False
    
    x = bracket(dt/bq, -bound, bound)
    kfun = kepler(x)
    
    past = dt < 0
    future = dt > 0
    
    upper = zeros_like(dt, dtype='float64')
    lower = zeros_like(dt, dtype='float64')
    oldx = zeros_like(dt, dtype='float64')
    
    lower[past] = x[past]
    upper[future] = x[future]
    
    while (kfun[past] > dt[past]).any():
        upper[past] = lower[past]
        lower[past] *= 2
        oldx[past] = x[past]
        x[past] = bracket(lower[past], -bound, bound)
        if (x[past] == oldx[past]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1}'
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))
        kfun[past] = kepler(x[past])
        
    while (kfun[future] < dt[future]).any():
        lower[future] = upper[future]
        upper[future] *= 2
        oldx[future] = x[future]
        x[future] = bracket(upper[future], -bound, bound)
        if (x[future] == oldx[future]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1} '
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))    
        kfun[future] = kepler(x[future])
    
    x = amin(array([upper, amax(array([lower, (lower+upper)/2]), axis=0)]), axis=0)
    
    lcount = zeros_like(dt)
    mostc = ones_like(dt)*1000
    not_done = (lower < x) * (x < upper)
    
    while not_done.any():
        kfun[not_done] = kepler(x[not_done])
        
        high = (kfun > dt) * not_done
        low = (kfun < dt) * not_done
        same = (~high * ~low) * not_done
        
        upper[high] = x[high]
        lower[low] = x[low]
        upper[same] = lower[same] = x[same]        
        
        condition = not_done * (mostc > 64) * (upper != 0) * (lower != 0)
        mostc[condition] = 64
        lcount[condition] = 0
        
        # vectorized version of min(upper, max(lower, (upper + lower)/2))        
        x[not_done] = amin(array([upper[not_done], amax(array([lower[not_done], (lower[not_done]+upper[not_done])/2]), axis=0)]), axis=0)
        lcount += 1
        not_done = (lower < x) * (x < upper) * (lcount < mostc)
    
    c0, c1, c2, c3 = stumpff(f*x*x)
    br = br0*c0 + x*(b2rv*c1 + x*(bq*c2))
    
    pc = 1 - qovr0 * x**2 * c2
    vc = dt - bq * x**3 * c3
    pcdot = -qovr0 / br * x * c1
    vcdot = 1 - bq / br * x**2 * c2
    
    if return_1d_array:
        position_prop = pc*position + vc*velocity
        velocity_prop = pcdot*position + vcdot*velocity
    else:
        position_prop = pc*tile(position[newaxis].T, dt.size) + vc*tile(velocity[newaxis].T, dt.size)
        velocity_prop = pcdot*tile(position[newaxis].T, dt.size) + vcdot*tile(velocity[newaxis].T, dt.size)
    
    return position_prop, velocity_prop
