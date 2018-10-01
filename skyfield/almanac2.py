from .api import load, Topos, EarthSatellite
from .constants import tau
from .optimizelib import newton, brent_min
from numpy import array, degrees, arcsin, where, diff, sort, hstack, linspace, ceil
from functools import partial
from scipy.misc import derivative

__all__ = ['meridian_transits', 'culminations', 'risings_settings', 
           'twilights', 'equinoxes', 'solstices', 'moon_quarters']

ts = load.timescale()

#%% Test Functions   
def _is_earth_based(location):
    """Returns True if ``location`` is geocentric or topocentric, False otherwise.
    """
    if not hasattr(location, 'positives'):
        return False
    elif location.target == 399:
        return True
    elif isinstance(location.positives[-1], Topos):
        return True
    else:
        return False
    
    
def _is_satellite(body):
    """Returns true if body or body.positives[-1] is an Earth Satellite
    """
    if isinstance(body, EarthSatellite):
        return True
    elif getattr(body, 'positives', ()) and isinstance(body.positives[-1], EarthSatellite):
        return True
    else:
        return False
    
    
#%% _find functions
def _find_value(f, value, partition_edges, slope_at_zero='positive', tol=1e-10):
    g = lambda t: (f(t) - value + 180) % 360 - 180
    g_edges = g(partition_edges)
            
    sign_changes = diff((g_edges>0).astype(int))
    if slope_at_zero=='positive':
        indices = where(sign_changes == 1)[0]
    elif slope_at_zero=='negative':
        indices = where(sign_changes == -1)[0]
    elif slope_at_zero=='any':
        indices = where(sign_changes != 0)[0]
        
    left_edges = partition_edges[indices]
    right_edges = partition_edges[indices+1]
            
    return newton(g, array(left_edges), array(right_edges), tol=tol)

def _find_extremes(f, partition_edges, find='min', tol=1e-15):
    if find == 'min':
        g = f
    elif find == 'max':
        g = lambda t:-f(t) + 360

    g_dot_edges = derivative(g, partition_edges, dx=1e-6)
            
    sign_changes = diff((g_dot_edges>0).astype(int))
    indices = where(sign_changes == 1)[0]
    left_edges = partition_edges[indices]
    right_edges = partition_edges[indices+1]
            
    return brent_min(g, array(left_edges), array(right_edges), tol=tol)


#%% Partial Functions

def _ra(observer, body, t):
    """Returns:
    observer.at(ts.tt(jd=t)).observe(body).apparent().radec(epoch='date')[0]._degrees
        
    or the right ascension of 'body' in degrees at terrestrial time 't' 
    when seen from 'observer'.
    """
    return observer.at(ts.tt(jd=t)).observe(body).apparent().radec(epoch='date')[0]._degrees


def _ecliptic_lon(observer, body, t):
    """Returns:
    observer.at(ts.tt(jd=t)).observe(body).apparent().ecliptic_latlon()[1].degrees
    or the ecliptic latitude of body in degrees at terrestrial time t.
    """
    if _is_earth_based(observer):
        return observer.at(ts.tt(jd=t)).observe(body).apparent().ecliptic_latlon(epoch='date')[1].degrees
    else:
        return observer.at(ts.tt(jd=t)).observe(body).ecliptic_latlon(epoch='date')[1].degrees


def _ecliptic_lon_diff(observer, body1, body2, t):
    """Returns the ecliptic longitudes of body1 minus that of body2
    in degrees at terrestrial time t. Ecliptic longitude is found with:
    observer.at(ts.tt(jd=t)).observe(body).apparent().ecliptic_latlon()[1].degrees
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

def meridian_transits(observer, body, t0, t1, kind='upper'):
    """Calculates data about upper and lower transits of the local meridian.
    
    This function searches between ``t0`` and ``t1`` for times when ``body``'s 
    local hour angle is 0h for an upper transit or 12h for a lower transit.

    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> sun = planets['sun']
    >>> earth = planets['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> transits(greenwich, mars, t0, t1)
    
    Arguments
    ---------
    observer : VectorSum
        VectorSum of earth + Topos
    body : Segment or VectorSum
        Vector representing the object whose transits are being found
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'upper'`` for upper transits
        ``'lower'`` for lower transits
        ``'all'`` for upper and lower transits
        
    Returns
    -------
    times : Time
        Times of transits
    """ 
    f = partial(_lha, observer, body)
    
    start = t0.tt
    end = t1.tt
    partition_width = .45
    num_partitions = int(ceil((end - start)/partition_width))
    partition_edges = linspace(start, end, num_partitions)
    
    if kind == 'all':
        upper_times = _find_value(f, 0, partition_edges)        
        lower_times = _find_value(f, 180, partition_edges)
        times = sort(hstack([upper_times, lower_times]))
    elif kind == 'upper':
        times = _find_value(f, 0, partition_edges)
    elif kind == 'lower':
        times = _find_value(f, 180, partition_edges)
    else:
        raise ValueError("kind must be 'all', 'upper', or 'lower'")
    
    return ts.tt(jd=times)


def culminations(observer, body, t0, t1, kind='upper'):
    """Calculates data about upper and lower culminations.
    
    This function searches between ``t0`` and ``t1`` for times when `body`'s 
    altitude reaches a local maximum or minimum. Finds upper culminations only 
    by default, but this can be changed with the ``kind`` keyword.

    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> sun = planets['sun']
    >>> earth = planets['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> culminations(greenwich, sun, t0, t1)
    
    Arguments
    ---------
    observer : VectorSum or Topos
        VectorSum of earth + Topos, or a plain Topos if ``body`` is a satellite
    body : Segment, VectorSum, or EarthSatellite
        Vector representing the object whose culminations are being found
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'upper'`` for local altitude maximums
        ``'lower'`` for local altitude minimums
        ``'all'`` for both upper and lower culminations
        
    Returns
    -------
    times : Time
        Times of altitude maximums or minimums
    """
    if _is_satellite(body):
        if getattr(body, 'positives', ()):
            body = body.positives[-1]
        if not isinstance(observer, Topos):
            observer = observer.positives[-1]
        f = partial(_satellite_alt, observer, body)
        period = tau/body.model.no/60/24 # days/orbit
        partition_width = period * .25
    else:
        f = partial(_alt, observer, body)
        partition_width = .25
    
    start = t0.tt
    end = t1.tt
    num_partitions = int(ceil((end - start)/partition_width))
    partition_edges = linspace(start, end, num_partitions)
    
    if kind == 'all':
        upper_times = _find_extremes(f, partition_edges, find='max')
        lower_times = _find_extremes(f, partition_edges)
        times = sort(hstack([upper_times, lower_times]))
    elif kind == 'upper':
        times = _find_extremes(f, partition_edges, find='max')
    elif kind == 'lower':
        times = _find_extremes(f, partition_edges)
    else:
        raise ValueError("kind must be 'all', 'upper', or 'lower'")

    return ts.tt(jd=times)
    

def risings_settings(observer, body, t0, t1, kind='all'):
    """Calculates data about when an object rises and sets.
    
    This function searches between ``t0`` and ``t1`` for times when `body`'s 
    altitude (uncorrected for refraction) is -34 arcminutes. The sun's and 
    moon's upper limb is used rather than their center.

    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> sun = planets['sun']
    >>> earth = planets['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> risings_settings(greenwich, sun, t0, t1, 'rise')
    
    Arguments
    ---------
    observer : VectorSum or Topos
        VectorSum of earth + Topos, or a plain Topos if ``body`` is a satellite
    body : Segment or VectorSum
        Vector representing the object whose rise/set times are being found
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'all'`` for when ``body`` rises and sets
        ``'rise'`` for when ``body`` rises
        ``'set'`` for when ``body`` sets
        
    Returns
    -------
    times : Time
        Times that `body` rises or sets
    """
    if _is_satellite(body):
        if getattr(body, 'positives', ()):
            body = body.positives[-1]
        if not isinstance(observer, Topos):
            observer = observer.positives[-1]
        f = partial(_satellite_alt, observer, body)
        value = -34/60
    elif body.target == 10: # sun
        f = partial(_alt, observer, body)
        value = -50/60
    elif body.target == 301: # moon
        f = partial(_moon_ul_alt, observer, body)
        value = -34/60
    else:
        f = partial(_alt, observer, body)
        value = 0
        
    
    
    start = t0.tt
    end = t1.tt
    body_culminations = culminations(observer, body, t0, t1, kind='all').tt
    partition_edges = hstack([start, body_culminations, end])    
    
    tol = 1e-15
    
    if kind == 'all':
        times = _find_value(f, value, partition_edges, slope_at_zero='any', tol=tol)
    elif kind == 'rise':
        times = _find_value(f, value, partition_edges, tol=tol)
    elif kind == 'set':
        times = _find_value(f, value, partition_edges, slope_at_zero='negative', tol=tol)
    else:
        raise ValueError("kind must be 'all', 'rise', or 'set'")
    
    return ts.tt(jd=times)
    

def twilights(observer, sun, t0, t1, kind='civil', begin_or_end='all'):
    """Calculates times when twilight starts in the morning.
   
    This function searches between ``t0`` and ``t1`` for times when the sun's 
    altitude is -6, -12, or -18 degrees for civil, nautical, or astronomical 
    twilights, respectively.

    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> earth = planets['earth']
    >>> greenwich = earth + Topos('51.5 N', '0 W')
    >>> t0 = ts.utc(2017, 1, 1)
    >>> t1 = ts.utc(2017, 1, 8)
    >>> twilight_ends(greenwich, t0, t1, kind='nautical')
    
    Arguments
    ---------
    observer : Segment or VectorSum
        Vector representing earth or a VectorSum of earth + Topos
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
    begin_or_end : str
        ``'all'`` for when twilight begins and ends
        ``'begin'`` for when twilight begins in the morning
        ``'end'`` for when twilight ends in the evening
        
    Returns
    -------
    times : Time
        Times that twilight starts in the morning and/ or ends in the evening
        """
        
    start = t0.tt
    end = t1.tt       
    sun_culminations = culminations(observer, sun, t0, t1, kind='all').tt
    partition_edges = hstack([start, sun_culminations, end])
    
    f = partial(_alt, observer, sun)
    if kind == 'civil':
        value = -6
    elif kind == 'nautical':
        value = -12
    elif kind == 'astronomical':
        value = -18
    else:
        raise ValueError("kind must be 'civil', 'nautical', or 'astronomical'")

    if begin_or_end == 'all':
        times = _find_value(f, value, partition_edges, slope_at_zero='any')
    elif begin_or_end == 'begin':
        times = _find_value(f, value, partition_edges)
    elif begin_or_end == 'end':
        times = _find_value(f, value, partition_edges, slope_at_zero='negative')
    else:
        raise ValueError("begin_or_end must be 'all', 'begin', or 'end'")

    return ts.tt(jd=times)
    
#%% Geocentric Phenomena, pg. 478 in Explanatory Supplement (1992)

def equinoxes(earth, t0, t1, kind='all'):
    """Calculates data about March and September equinoxes.
    
    This function searches between ``t0`` and ``t1`` for times when the sun's 
    ecliptic longitude  is 0 and 180 degrees, for vernal and autumnal 
    equinoxes, respectively.
    
    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> earth = planets['earth']
    >>> t0 = ts.utc(2017)
    >>> t1 = ts.utc(2018)
    >>> equinoxes(earth, t0, t1, kind='march')
    
    Arguments
    ---------
    earth : Segment
        Vector representing earth
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'march'``  finds equinoxes in march
        ``'september'``  finds equinoxes in september
        ``'all'``  finds all equinoxes
        
    Returns
    -------
    times : Time
        Times of solstices
    """
    sun = earth.ephemeris['sun']

    f = partial(_ecliptic_lon, earth, sun)

    start = t0.tt
    end = t1.tt
    partition_width = 365*.45
    num_partitions = int(ceil((end - start)/partition_width))
    partition_edges = linspace(start, end, num_partitions)
    
    if kind == 'all':
        march_times = _find_value(f, 0, partition_edges)
        september_times = _find_value(f, 180, partition_edges)
        times = sort(hstack([march_times, september_times]))
    elif kind == 'march':
        times = _find_value(f, 0, partition_edges)
    elif kind == 'september':
        times = _find_value(f, 180, partition_edges)
    else:
        raise ValueError("kind must be 'all', 'march', or 'september'")
        
    return ts.tt(jd=times)
        

def solstices(earth, t0, t1, kind='all'):
    """Calculates data about June and December solstices.
    
    This function searches between ``t0`` and ``t1`` for times when the sun's 
    ecliptic longitude is 90 and 270 degrees, for June and December solstices, 
    respectively.
    
    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> earth = planets['earth']
    >>> t0 = ts.utc(2017, 1)
    >>> t1 = ts.utc(2017, 2)
    >>> solstices(earth, t0, t1, kind='June')
    
    Arguments
    ---------
    earth : Segment
        Vector representing earth
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'june'``  finds solstices in June
        ``'december'``  finds solstices in December
        ``'all'``  finds all solstices
        
    Returns
    -------
    times : Time
        Times of solstices
    """
    sun = earth.ephemeris['sun']
    
    f = partial(_ecliptic_lon, earth, sun)
    
    start = t0.tt
    end = t1.tt
    partition_width = 365*.45
    num_partitions = int(ceil((end - start)/partition_width))
    partition_edges = linspace(start, end, num_partitions)
    
    if kind == 'all':
        june_times = _find_value(f, 90, partition_edges)
        december_times = _find_value(f, 270, partition_edges)
        times = sort(hstack([june_times, december_times]))
    elif kind == 'june':
        times = _find_value(f, 90, partition_edges)
    elif kind == 'december':
        times = _find_value(f, 270, partition_edges)
    else:
        raise ValueError("kind must be 'all', 'june', or 'december'")
    
    return ts.tt(jd=times)


def moon_quarters(moon, t0, t1, kind='all'):
    """Calculates data about new and full moons, first and last quarters.
    
    This function searches between ``t0`` and ``t1`` for times when the moon's 
    geocentric ecliptic longitude minus that of the sun is:
        
        * 0 degrees for new moon
        * 90 degrees for first quarter
        * 180 degrees for full moon
        * 270 degrees for last quarter 
    
    Example
    -------
    >>> planets = load('de430t.bsp')
    >>> moon = planets['moon']
    >>> t0 = ts.utc(2017, 1)
    >>> t1 = ts.utc(2017, 2)
    >>> moon_quarters(moon, t0, t1, kind='full')
    
    Arguments
    ---------
    moon : Segment
        Vector representing the moon
    t0 : Time
        Time object of length 1 representing the start of the search interval
    t1 : Time
        Time object of length 1 representing the end of the search interval
    kind : str
        ``'new'``  for new moons
        ``'first'``  for first quarters
        ``'full'``  for full moons
        ``'last'``  for last quarters
        ``'all'``  for all quarters of the moon
        
    Returns
    -------
    times : Time
        Times of moon quarters
    """ 
    earth = moon.ephemeris['earth']
    sun = moon.ephemeris['sun']
    
    f = partial(_ecliptic_lon_diff, earth, moon, sun)
    
    start = t0.tt
    end = t1.tt    
    partition_width = 29*.25
    num_partitions = int(ceil((end - start)/partition_width))
    partition_edges = linspace(start, end, num_partitions)
    
    if kind == 'all':
        new_times = _find_value(f, 0, partition_edges)
        first_times = _find_value(f, 90, partition_edges)
        full_times = _find_value(f, 180, partition_edges)
        last_times = _find_value(f, 270, partition_edges)
        times = sort(hstack([new_times, first_times, full_times, last_times]))
    elif kind == 'new':
        times = _find_value(f, 0, partition_edges)
    elif kind == 'first':
        times = _find_value(f, 90, partition_edges)
    elif kind == 'full':
        times = _find_value(f, 180, partition_edges)
    elif kind == 'last':
        times = _find_value(f, 270, partition_edges)
    else:
        raise ValueError("kind must be 'all', 'new', 'first', 'full', or 'last'")

    return ts.tt(jd=times)
