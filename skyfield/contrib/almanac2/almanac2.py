"""
Rough outline of how this module finds times of phenomena:

    * Split the time interval evenly into partitions. The partitions 
    should be small enough that no partitions contain multiple target values.
    * Evaluate the objective function at the partition edges to determine which 
    partitions contain the points of interest.
    * Pass the partitions that contain points of interest to the vectorized 
    optimization routines.
    
Some other things to be aware of are:
    
    * The variable ``f`` is the plain objective function that returns numbers 
    directly from skyfield.
    * The variable ``g`` is ``f`` transformed such that either all the target 
    values appear to the secant method as roots, or such that all of the 
    extremes appear to Brent's method as minima.
    
"""
from skyfield.api import load, Topos, EarthSatellite, Angle, Star
from skyfield.constants import tau
from skyfield.vectorlib import VectorSum
from optimizelib import secant, brent_min
from numpy import (degrees, arcsin, isfinite, hstack, nan, empty, linspace, 
                   ceil, sign, nonzero, zeros, empty_like, array)
from functools import partial

__all__ = ['meridian_transits', 'culminations', 'risings_settings', 
           'twilights', 'seasons', 'moon_phases']

ts = load.timescale()
    
#%% functions for finding/ isolating partitions
            
def _find_value(f, values, partition_edges, slope='any'):
    """ Evaluates ``f`` at ``partition_edges`` and returns the left and right 
    edges of those partitions that contain ``values``.
    
    Arguments
    ---------
    f : function
        Objective function. Must accept and return ndarrays.
    values : list
        List of values of ``f`` that are to be found
    partition_edges : ndarray
        Array of partition edges
    slope : str
        'positive' to find ``values`` when the slope is positive
        'negative' to find ``values`` when the slope is negative
        'any' to find ``values`` at any slope
        
    Returns
    -------
    jd0 : ndarray
        Left edges of those partitions found to contain ``values``
    jd1 : ndarray
        Right edges of those partitions found to contain ``values``
    targets : ndarray
        target from ``values`` that each partition is found to contain
    f0 : ndarray
        ``f(jd0)``
    f1 : ndarray
        ``f(jd1)``
    is_positive : ndarray, dtype=bool
        whether or not the slope of ``f`` in the partition is positive
    """
    f_values = f(partition_edges)             
                
    targets = empty(len(partition_edges)-1)
    targets.fill(nan)
    
    for value in values:
        g_values = (f_values - value + 180)%360 - 180
        if slope=='positive':
            found = (sign(g_values[:-1])==-1) * (sign(g_values[1:])==1)
        elif slope=='negative':
            found = (sign(g_values[:-1])==1) * (sign(g_values[1:])==-1)
        elif slope=='any':
            found = sign(g_values[:-1]) != sign(g_values[1:])
        
        if (isfinite(targets) * found).any():
            raise ValueError('Multiple target values found in the same partition. Make the partitions smaller.')
            
        targets[found] = value
        
    indices = nonzero(isfinite(targets))[0]
    jd0 = partition_edges[indices]
    jd1 = partition_edges[indices+1]
    f0 = f_values[indices]
    f1 = f_values[indices+1]
    
    is_positive = ((sign(g_values[:-1])==-1) * (sign(g_values[1:])==1))[indices]
    
    return jd0, jd1, targets[indices], f0, f1, is_positive


def _find_extremes(f, partition_edges, find='min'):
    """ Evaluates the derivative of ``f`` at ``partition_edges`` and returns 
    the left and right edges of the partitions that contain extrema.
    
    Arguments
    ---------
    f : function
        Objective function. Must accept and return ndarrays.
    partition_edges : ndarray
        Array of partition edges
    find : str
        'min' to find minima
        'max' to find maxima
        'any' to find all extrema
        
    Returns
    -------
    jd0 : ndarray
        Left edges of those partitions found to contain ``values``
    jd1 : ndarray
        Right edges of those partitions found to contain ``values``
    targets : ndarray
        target from ``values`` that each partition is found to contain
    f0 : ndarray
        ``f(jd0)``
    f1 : ndarray
        ``f(jd1)``
    minimum : ndarray, dtype=bool
        whether or not each partition contains a minima (False means the 
        partition contains a maxima)
    """  
    # evaluate the derivative using forward difference method. Combine both 
    # sides of the approximation so only one function call is necessary.
    step_size = 1e-6
    combined_array = hstack([partition_edges, partition_edges+step_size])
    combined_results = f(combined_array)
    left = combined_results[:len(partition_edges)]
    right = combined_results[len(partition_edges):]
    f_dot = (right - left)/step_size
    
    has_extreme = zeros(len(partition_edges)-1)
    
    if find == 'min' or find == 'any':
        has_min = (sign(f_dot[:-1])==-1) * (sign(f_dot[1:])==1)
        has_extreme[has_min] = -1
    
    if find == 'max' or find == 'any':
        has_max = (sign(f_dot[:-1])==1) * (sign(f_dot[1:])==-1)
        if find == 'any' and (has_min * has_max).any():
            raise ValueError('Multiple extremes found in the same partition. Make the partitions smaller.')
        has_extreme[has_max] = 1
        
    indices = nonzero(has_extreme)[0]
    jd0 = partition_edges[indices]
    jd1 = partition_edges[indices+1]
    minimum = has_extreme[indices]==-1
    f0 = left[indices]
    f1 = left[indices+1]
    
    return jd0, jd1, minimum, f0, f1


def divide_evenly(start, end, max_width):
    """ Evenly divides the interval between start and end into intervals that 
    are at most max_width wide.
    
    Arguments
    ---------
    start : float
        Start of the interval
    end : float
        End of the interval
    max_width : float
        Maximum width of the divisions
        
    Returns
    -------
    divisions : ndarray
        Resulting array
    """
    num_partitions = int(ceil((end - start)/max_width))
    return linspace(start, end, num_partitions+1)


#%% Objective Functions

def _ra(observer, body, t):
    """Returns the right ascension of 'body' in degrees at terrestrial time 't' 
    when seen from 'observer'.
    """
    return observer.at(ts.tt(jd=t)).observe(body).apparent().radec(epoch='date')[0]._degrees


def _ecliptic_lon(observer, body, t):
    """Returns the ecliptic latitude of body in degrees at terrestrial time t.
    """
    return observer.at(ts.tt(jd=t)).observe(body).apparent().ecliptic_latlon(epoch='date')[1].degrees


def _ecliptic_lon_diff(observer, body1, body2, t):
    """Returns the ecliptic longitudes of body1 minus that of body2
    in degrees at terrestrial time t.
    """
    diff = _ecliptic_lon(observer, body1, t) - _ecliptic_lon(observer, body2, t)
    return (diff + 180)%360 - 180


def _local_sidereal(observer, t):
    """Returns observer's local apparent sidereal time at t in degrees"""
    return ts.tt(jd=t).gast*15 + observer.positives[-1].longitude.degrees


def _lha(observer, body, t):
    """Returns the local hour angle of `body` in degrees when seen from 
    `observer` at terrestrial time `t`.
    """
    return _local_sidereal(observer, t) - _ra(observer, body, t)


def _alt(observer, body, t):
    """Returns the altitude of `body` in degrees when seen from `observer` at 
    terrestrial time `t`.
    """
    return observer.at(ts.tt(jd=t)).observe(body).apparent().altaz()[0].degrees


def _moon_ul_alt(observer, moon, t):
    """Returns the altitude of the moon's upper limb in degrees when seen from
    observer at terrestrial time t.
    
    Calculates the moon's apparent semidiameter from the location of 
    ``observer``
    """
    moon_distance = observer.at(ts.tt(jd=t)).observe(moon).apparent().distance().au
    moon_radius = 1.161781e-5 #au
    moon_sd = degrees(arcsin(moon_radius/moon_distance))
    return _alt(observer, moon, t) + moon_sd


def _satellite_alt(observer, satellite, t):
    """Returns the altitude of `satellite` in degrees when seen from `observer` 
    at terrestrial time `t`.
    """
    return (satellite - observer).at(ts.tt(jd=t)).altaz()[0].degrees

#%% Topocentric Phenomena, pg. 482 in Explanatory Supplement (1992)

def meridian_transits(observer, body, t0, t1):
    """Calculates times of upper and lower transits of the local meridian.
    
    This function searches between ``t0`` and ``t1`` for times when ``body``'s 
    local hour angle is 0h for an upper transit or 12h for a lower transit.

    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> sun = ephem['sun']
    >>> earth = ephem['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> times, hour_angles = transits(greenwich, mars, t0, t1)
    >>>
    >>> upper_transits = times[hour_angles.hours==0]
    >>> lower_transits = times[hour_angles.hours==12]
    
    Arguments
    ---------
    observer : VectorSum
        VectorSum of earth + Topos
    body : Segment or VectorSum
        Vector representing the object whose transits are being found.
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
        
    Returns
    -------
    times : Time
        Times of transits
    hour_angles : Angle
        Local Hour Angle of ``body`` at ``times``
    """ 
    if isinstance(observer, Topos):
        raise ValueError('`observer` should be a VectorSum of earth + Topos.')
    
    if isinstance(body, EarthSatellite) or (isinstance(body, VectorSum) and isinstance(body.positives[-1], EarthSatellite)):
        raise ValueError("meridian_transits doesn't support EarthSatellites.")
    
    f = partial(_lha, observer, body)    
    partition_edges = divide_evenly(t0.tt, t1.tt, .2)
    
    jd0, jd1, targets, f0, f1, _ = _find_value(f, [0, 180], partition_edges, slope='positive')
    
    times = secant(f, jd0, jd1, targets, f0, f1)
    
    return ts.tt(jd=times), Angle(degrees=targets, preference='hours')


def culminations(observer, body, t0, t1):
    """Calculates times of upper and lower culminations.
    
    This function searches between ``t0`` and ``t1`` for times when `body`'s 
    altitude reaches a local maximum or minimum. Finds upper culminations only 
    by default, but this can be changed with the ``kind`` keyword.

    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> sun = ephem['sun']
    >>> earth = planets['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> times, kinds = culminations(greenwich, sun, t0, t1)
    >>>
    >>> upper_culminations = times[kinds=='upper']
    >>> lower_culminations = times[kinds=='lower']
    
    Arguments
    ---------
    observer : VectorSum or Topos
        VectorSum of earth + Topos. If ``body`` is an EarthSatellite, 
        ``observer`` can be either a VectorSum or a plain Topos object.
    body : Segment, VectorSum, or EarthSatellite
        Vector representing the object whose culminations are being found. For 
        EarthSatellites use a plain Earthsatellite and not a VectorSum.
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
        
    Returns
    -------
    times : Time
        Times of altitude maximums or minimums
    kinds : ndarray, dtype=str
        array containing 'upper' for upper culminations, or 'lower' for lower 
        culminations
    """    
    if isinstance(body, EarthSatellite) and isinstance(observer, VectorSum):
        observer = observer.positives[-1]
        
    if not isinstance(body, EarthSatellite) and not isinstance(observer, VectorSum):
        raise ValueError('Unless ``body`` is an EarthSatellite, ``observer``'
                         'must be a VectorSum of Earth + Topos')
        
    if isinstance(body, VectorSum) and isinstance(body.positives[-1], EarthSatellite):
        raise ValueError('`body` should be a plain EarthSatellite, not a VectorSum')
    
    if isinstance(body, EarthSatellite):
        f = partial(_satellite_alt, observer, body)
        period = tau/body.model.no/60/24 # days/orbit
        max_width = period * .2
    else:
        f = partial(_alt, observer, body)
        max_width = .2
    
    partition_edges = divide_evenly(t0.tt, t1.tt, max_width)

    jd0, jd1, minimum, f0, f1 = _find_extremes(f, partition_edges, 'any')

    if len(jd0)!=0:
        times = brent_min(f, jd0, jd1, minimum, f0, f1, tol=1e-15)
    else:
        times = array([])

    kinds = empty_like(minimum, dtype='U5')
    kinds[minimum] = 'lower'
    kinds[~minimum] = 'upper'

    return ts.tt(jd=times), kinds
    

def risings_settings(observer, body, t0, t1):
    """Calculates times when an object rises and sets.
    
    This function searches between ``t0`` and ``t1`` for times when `body`'s 
    altitude (uncorrected for refraction) is -34 arcminutes. The sun's and 
    moon's upper limb is used rather than their center.

    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> sun = ephem['sun']
    >>> earth = ephem['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> times, kinds = risings_settings(greenwich, sun, t0, t1)
    >>>
    >>> risings = times[kinds=='rise']
    >>> settings = times[kinds=='set']
    
    Arguments
    ---------
    observer : VectorSum or Topos
        VectorSum of earth + Topos. If ``body`` is an EarthSatellite, 
        ``observer`` can be either a VectorSum or a plain Topos object.
    body : Segment, VectorSum, or EarthSatellite
        Vector representing the object whose rise/set times are being found. 
        For EarthSatellite objects use a plain EarthSatellite and not a 
        VectorSum.
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
        
    Returns
    -------
    times : Time
        Times that `body` rises or sets
    kinds : ndarray, dtype=str
        array containing 'rise' for risings, or 'set' for settings
    """
    if isinstance(body, EarthSatellite) and isinstance(observer, VectorSum):
        observer = observer.positives[-1]
    
    if not isinstance(body, EarthSatellite) and not isinstance(observer, VectorSum):
        raise ValueError('Unless ``body`` is an EarthSatellite, ``observer``'
                         'must be a VectorSum of Earth + Topos')
        
    if isinstance(body, VectorSum) and isinstance(body.positives[-1], EarthSatellite):
        raise ValueError('`body` should be a plain EarthSatellite, not a VectorSum')
    
    if isinstance(body, EarthSatellite):
        f = partial(_satellite_alt, observer, body)
        value = [-34/60]
    elif isinstance(body, Star):
        f = partial(_alt, observer, body)
        value = [0]
    elif body.target == 10: # sun
        f = partial(_alt, observer, body)
        value = [-50/60]
    elif body.target == 301: # moon
        f = partial(_moon_ul_alt, observer, body)
        value = [-34/60]
    else:
        f = partial(_alt, observer, body)
        value = [0]
    
    body_culminations = culminations(observer, body, t0, t1)[0].tt
    partition_edges = hstack([t0.tt, body_culminations, t1.tt])    
    
    jd0, jd1, targets, f0, f1, is_positive = _find_value(f, value, partition_edges)
    
    times = secant(f, jd0, jd1, targets, f0, f1, tol=1e-15)
    
    kinds = empty_like(is_positive, dtype='U4')
    kinds[is_positive] = 'rise'
    kinds[~is_positive] = 'set'

    return ts.tt(jd=times), kinds
    

def twilights(observer, sun, t0, t1, kind='civil'):
    """Calculates times when twilight starts or ends in the morning or evening.
   
    This function searches between ``t0`` and ``t1`` for times when the sun's 
    altitude is -6, -12, or -18 degrees for civil, nautical, or astronomical 
    twilights, respectively.

    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> earth = ephem['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> times, am_pm = (greenwich, t0, t1, kind='civil')
    >>>
    >>> am_twilights = times[am_pm=='am']
    >>> pm_twilights = times[am_pm=='pm']
    
    Arguments
    ---------
    observer : VectorSum
        VectorSum of earth + Topos
    sun : Segment
        Segment representing the sun
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'civil'`` for civil twilight
        ``'nautical'`` for nautical twilight
        ``'astronomical'`` for astronomical twilight
        
    Returns
    -------
    times : Time
        Times that twilight starts in the morning and/ or ends in the evening
    am_pm : ndarray, dtype=str
        array containing 'am' for morning twilight, or 'pm' for evening twilight
        """
    if not isinstance(observer.positives[-1], Topos):
        raise ValueError('`observer` should be a VectorSum of earth + Topos.')
        
    sun_culminations = culminations(observer, sun, t0, t1)[0].tt
    partition_edges = hstack([t0.tt, sun_culminations, t1.tt])
    
    f = partial(_alt, observer, sun)
    if kind == 'civil':
        value = [-6]
    elif kind == 'nautical':
        value = [-12]
    elif kind == 'astronomical':
        value = [-18]
    else:
        raise ValueError("kind must be 'civil', 'nautical', or 'astronomical'")

    jd0, jd1, targets, f0, f1, is_positive = _find_value(f, value, partition_edges)

    times = secant(f, jd0, jd1, targets, f0, f1)

    am_pm = empty_like(is_positive, dtype='U2')
    am_pm[is_positive] = 'am'
    am_pm[~is_positive] = 'pm'

    return ts.tt(jd=times), am_pm
    
#%% Geocentric Phenomena, pg. 478 in Explanatory Supplement (1992)

def seasons(earth, t0, t1):
    """Calculates times of equinoxes and solstices.
    
    This function searches between ``t0`` and ``t1`` for times when the sun's 
    ecliptic longitude is:
        
        * 0 degrees for march equinoxes
        * 90 degrees for june solstices
        * 180 degrees for september equinoxes
        * 270 degrees for december solstices
    
    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> earth = ephem['earth']
    >>> t0 = ts.utc(2017)
    >>> t1 = ts.utc(2018)
    >>> times, lons = seasons(earth, t0, t1)
    >>>
    >>> march_equinoxes = times[lons.degrees==0]
    >>> june_solstices = times[lons.degrees==90]
    >>> sept_equinoxes = times[lons.degrees==180]
    >>> dec_solstices = times[lons.degrees==270]
    
    Arguments
    ---------
    earth : Segment
        Vector representing earth
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
        
    Returns
    -------
    times : Time
        Times of solstices
    longitudes : Angle
        sun's ecliptic longitude at ``times``
    """
    sun = earth.ephemeris['sun']

    f = partial(_ecliptic_lon, earth, sun)

    partition_edges = divide_evenly(t0.tt, t1.tt, 365*.2)

    jd0, jd1, targets, f0, f1, _ = _find_value(f, [0, 90, 180, 270], partition_edges, slope='positive')
        
    times = secant(f, jd0, jd1, targets, f0, f1)

    return ts.tt(jd=times), Angle(degrees=targets)


def moon_phases(moon, t0, t1):
    """Calculates times of the phases of the moon.
    
    This function searches between ``t0`` and ``t1`` for times when the moon's 
    geocentric ecliptic longitude minus that of the sun is:
        
        * 0 degrees for new moon
        * 90 degrees for first quarter
        * 180 degrees for full moon
        * 270 degrees for last quarter 
    
    Example
    -------
    >>> ephem = load('de430t.bsp')
    >>> moon = ephem['moon']
    >>> t0 = ts.utc(2017, 1)
    >>> t1 = ts.utc(2017, 2)
    >>> times, lon_diffs = moon_phases(moon, t0, t1)
    >>>
    >>> new_moons = times[lon_diffs.degrees==0]
    >>> first_quarters = times[lon_diffs.degrees==90]
    >>> full_moons = times[lon_diffs.degrees==180]
    >>> last_quarters = times[lon_diffs.degrees==270]
    
    Arguments
    ---------
    moon : Segment
        Vector representing the moon
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
        
    Returns
    -------
    times : Time
        Times of moon quarters
    longitude_diffs: Angle
        moon's ecliptic longitude - sun's ecliptic longitude at ``times``
        
    """ 
    earth = moon.ephemeris['earth']
    sun = moon.ephemeris['sun']
    
    f = partial(_ecliptic_lon_diff, earth, moon, sun)
    
    partition_edges = divide_evenly(t0.tt, t1.tt, 29*.2)

    jd0, jd1, targets, f0, f1, _ = _find_value(f, [0, 90, 180, 270], partition_edges, slope='positive')
        
    times = secant(f, jd0, jd1, targets, f0, f1)

    return ts.tt(jd=times), Angle(degrees=targets)
