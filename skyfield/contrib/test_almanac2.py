from skyfield.api import Loader, Topos
from skyfield.timelib import Time
from almanac2 import (meridian_transits, culminations, twilights, 
                       risings_settings, seasons, moon_quarters)
from almanac2 import (_ecliptic_lon_diff, _moon_ul_alt, _lha, _alt,
                       _satellite_alt, _ecliptic_lon)
from numpy import concatenate, ndarray
from functools import partial

load = Loader(r'C:\Users\Josh\Scripts\Skyfield_Data')

ts = load.timescale()
ephem = load('de430.bsp')
earth = ephem['earth']
sun = ephem['sun']
moon = ephem['moon']
mars = ephem['mars barycenter']
venus = ephem['venus']
jupiter_bc = ephem['jupiter barycenter']
mercury = ephem['mercury']
greenwich = earth + Topos('51.5 N', '0 W')

satellites = load.tle(r'stations.txt')
ISS = earth + satellites['ISS']

jup_ephem = load('jup310.bsp')
jupiter = jup_ephem['jupiter']
io = jup_ephem['io']

minute = 1/24/60
sec = minute/60
ms = sec/1000


def compare(value, expected_value, epsilon):
    if hasattr(value, '__len__') or hasattr(expected_value, '__len__'):
        assert max(abs(value - expected_value)) <= epsilon
    else:
        assert abs(value - expected_value) <= epsilon


def is_root(f, times, target_f, epsilon):
    left_times = times - epsilon
    right_times = times + epsilon
    left_f = (f(left_times) - target_f + 180)%360 - 180
    right_f = (f(right_times) - target_f + 180)%360 - 180
    
    # TODO: could this be success = sign(left_edge) != sign(right_edge)
    success_increasing = (left_f <= 0) * (right_f >= 0)
    success_decreasing = (left_f >= 0) * (right_f <= 0)
    success = success_increasing + success_decreasing
    
    if isinstance(success, ndarray):
        return success.all()
    else:
        return success
    
    
def is_extreme(f, times, epsilon):
    left_times = times - epsilon
    right_times = times + epsilon
    center_f = f(times)
    left_f = (f(left_times) - center_f + 180)%360 - 180
    right_f = (f(right_times) - center_f + 180)%360 - 180
    
    # TODO: could this be success = sign(left_edge) == sign(right_edge)
    is_max = (left_f <= 0) * (right_f <= 0)
    is_min = (left_f >= 0) * (right_f >= 0)
    success = is_max + is_min
    
    if isinstance(success, ndarray):
        return success.all()
    else:
        return success


# Data Source:
# web.archive.org/web/20171223014917/http://aa.usno.navy.mil/data/docs/EarthSeasons.php
def test_seasons():
    t0 = ts.utc(2000)
    t1 = ts.utc(2026)
    times = seasons(earth, t0, t1)
    
    assert isinstance(times, Time)
    assert len(times) == 104
    
    compare(times[0].tt, ts.utc(2000, 3, 20, 7, 35).tt, minute/2)
    compare(times[1].tt, ts.utc(2000, 6, 21, 1, 48).tt, minute/2)
    
    # Check that the parts equal the whole
    march_times = seasons(earth, t0, t1, 'march')
    june_times = seasons(earth, t0, t1, 'june')
    sept_times = seasons(earth, t0, t1, 'september')
    dec_times = seasons(earth, t0, t1, 'december')
    all_times = concatenate([march_times.tt, june_times.tt, sept_times.tt, dec_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data
    f = partial(_ecliptic_lon, earth, sun)
    assert is_root(f, march_times.tt, 0, ms/2)
    assert is_root(f, june_times.tt, 90, ms/2)
    assert is_root(f, sept_times.tt, 180, ms/2)   
    assert is_root(f, dec_times.tt, 270, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_phases.pl?year=2017&month=1&day=1&nump=50&format=p
def test_moon_quarters():
    t0 = ts.utc(2017)
    t1 = ts.utc(2018)
    times = moon_quarters(moon, t0, t1)
    
    assert isinstance(times, Time)
    assert len(times) == 49
    
    compare(times[0].tt, ts.utc(2017, 1, 5, 19, 47).tt, minute/2)
    
    new_times = moon_quarters(moon, t0, t1, 'new')
    first_times = moon_quarters(moon, t0, t1, 'first')
    full_times = moon_quarters(moon, t0, t1, 'full')
    last_times = moon_quarters(moon, t0, t1, 'last')
    all_times = concatenate([new_times.tt, first_times.tt, full_times.tt, last_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data    
    f = partial(_ecliptic_lon_diff, earth, moon, sun)
    assert is_root(f, new_times.tt, 0, ms/2)
    assert is_root(f, first_times.tt, 90, ms/2)
    assert is_root(f, full_times.tt, 180, ms/2)
    assert is_root(f, last_times.tt, 270, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=0&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_sun_risings_settings():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = risings_settings(greenwich, sun, t0, t1)
    
    assert isinstance(times, Time)
    assert len(times) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 8, 6).tt, minute/2)
    
    rise_times = risings_settings(greenwich, sun, t0, t1, 'rise')
    set_times = risings_settings(greenwich, sun, t0, t1, 'set')
    all_times = concatenate([rise_times.tt, set_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data    
    f = partial(_alt, greenwich, sun)
    assert is_root(f, times.tt, -50/60, ms/2)


# Data Source:   
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=1&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_moon_risings_settings():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = risings_settings(greenwich, moon, t0, t1)
    
    assert isinstance(times, Time)
    assert len(times) == 60
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 9, 46).tt, minute/2)
    
    rise_times = risings_settings(greenwich, moon, t0, t1, 'rise')
    set_times = risings_settings(greenwich, moon, t0, t1, 'set')
    all_times = concatenate([rise_times.tt, set_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all() 
    
    # Check that the found times produce the correct data
    f = partial(_moon_ul_alt, greenwich, moon)
    assert is_root(f, times.tt, -34/60, ms/2)
    
    
#TODO: how can this test be stable as tle's change?
def test_ISS_risings_settings():
    t0 = ts.utc(2017, 6, 1)
    t1 = ts.utc(2017, 6, 2)
    times = risings_settings(greenwich, ISS, t0, t1)
    
    assert isinstance(times, Time)
    assert len(times) == 13
    
    rise_times = risings_settings(greenwich, ISS, t0, t1, 'rise')
    set_times = risings_settings(greenwich, ISS, t0, t1, 'set')
    all_times = concatenate([rise_times.tt, set_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all() 
    
    # Check that the found times produce the correct data
    f = partial(_satellite_alt, greenwich, ISS)
    assert is_root(f, times.tt, -34/60, ms/2)
    
    # Check that same result is found if plain Topos is used as observer
    times2 = risings_settings(greenwich.positives[-1], ISS, t0, t1)
    assert (times.tt == times2.tt).all()
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=2&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_civil_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = twilights(greenwich, sun, t0, t1, 'civil')

    assert isinstance(times, Time)
    assert len(times) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 7, 26).tt, minute/2)
    
    begin_times = twilights(greenwich, sun, t0, t1, 'civil', 'begin')
    end_times = twilights(greenwich, sun, t0, t1, 'civil', 'end')
    all_times = concatenate([begin_times.tt, end_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data
    f = partial(_alt, greenwich, sun)
    assert is_root(f, begin_times.tt, -6, ms/2)
    assert is_root(f, end_times.tt, -6, ms/2)


# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=3&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_nautical_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = twilights(greenwich, sun, t0, t1, 'nautical')
    
    assert isinstance(times, Time)
    assert len(times) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 6, 43).tt, minute/2)
    
    begin_times = twilights(greenwich, sun, t0, t1, 'nautical', 'begin')
    end_times = twilights(greenwich, sun, t0, t1, 'nautical', 'end')
    all_times = concatenate([begin_times.tt, end_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all() 
    
    # Check that the found times produce the correct data
    f = partial(_alt, greenwich, sun)
    assert is_root(f, begin_times.tt, -12, ms/2)
    assert is_root(f, end_times.tt, -12, ms/2)
    
    
# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=4&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_astronomical_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = twilights(greenwich, sun, t0, t1, 'astronomical')
    
    assert isinstance(times, Time)
    assert len(times) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 6, 2).tt, minute/2)
    
    begin_times = twilights(greenwich, sun, t0, t1, 'astronomical', 'begin')
    end_times = twilights(greenwich, sun, t0, t1, 'astronomical', 'end')
    all_times = concatenate([begin_times.tt, end_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all() 
    
    # Check that the found times produce the correct data
    f = partial(_alt, greenwich, sun)
    assert is_root(f, begin_times.tt, -18, ms/2)
    assert is_root(f, end_times.tt, -18, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_mrst2.pl?form=2&ID=AA&year=2017&month=1&day=1&reps=31&body=10&place=Greenwich&lon_sign=-1&lon_deg=&lon_min=&lon_sec=&lat_sign=1&lat_deg=51&lat_min=30&lat_sec=&height=&tz=&tz_sign=-1
def test_sun_transits():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = meridian_transits(greenwich, sun, t0, t1, 'all')
    
    assert isinstance(times, Time)
    assert len(times) == 62
        
    compare(times[1].tt, ts.utc(2017, 1, 1, 12, 4).tt, minute/2)
    
    upper_times = meridian_transits(greenwich, sun, t0, t1, 'upper')
    lower_times = meridian_transits(greenwich, sun, t0, t1, 'lower')
    all_times = concatenate([upper_times.tt, lower_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()

    # Check that the found times produce the correct data
    f = partial(_lha, greenwich, sun)
    assert is_root(f, upper_times.tt, 0, ms/2)
    assert is_root(f, lower_times.tt, 180, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_mrst2.pl?form=2&ID=AA&year=2017&month=1&day=1&reps=31&body=4&place=Greenwich&lon_sign=-1&lon_deg=&lon_min=&lon_sec=&lat_sign=1&lat_deg=51&lat_min=30&lat_sec=&height=&tz=&tz_sign=-1
def test_mars_transits():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = meridian_transits(greenwich, mars, t0, t1, 'all')
    
    assert isinstance(times, Time)
    assert len(times) == 62
    
    compare(times[1].tt, ts.utc(2017, 1, 1, 16, 2).tt, minute/2)
    
    upper_times = meridian_transits(greenwich, mars, t0, t1, 'upper')
    lower_times = meridian_transits(greenwich, mars, t0, t1, 'lower')
    all_times = concatenate([upper_times.tt, lower_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data
    f = partial(_lha, greenwich, mars)
    assert is_root(f, upper_times.tt, 0, ms/2)
    assert is_root(f, lower_times.tt, 180, ms/2)
    

def test_sun_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = culminations(greenwich, sun, t0, t1, 'all')
    
    assert isinstance(times, Time)
    assert len(times) == 62
    
    upper_times = culminations(greenwich, sun, t0, t1, 'upper')
    lower_times = culminations(greenwich, sun, t0, t1, 'lower')
    all_times = concatenate([upper_times.tt, lower_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()
    
    # Check that the found times produce the correct data
    f = partial(_alt, greenwich, sun)        
    assert is_extreme(f, times.tt, 2*ms)
    
    
def test_moon_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times = culminations(greenwich, moon, t0, t1, 'all')
    
    assert isinstance(times, Time)
    assert len(times) == 60
    
    upper_times = culminations(greenwich, moon, t0, t1, 'upper')
    lower_times = culminations(greenwich, moon, t0, t1, 'lower')
    all_times = concatenate([upper_times.tt, lower_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()

    # Check that the found times produce the correct data
    f = partial(_alt, greenwich, moon)    
    assert is_extreme(f, times.tt, 6*ms)
    
    
#TODO: how can this test be stable as tle's change?
def test_ISS_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 2)
    times = culminations(greenwich, ISS, t0, t1, 'all')
    
    assert isinstance(times, Time)
    assert len(times) == 31

    upper_times = culminations(greenwich, ISS, t0, t1, 'upper')
    lower_times = culminations(greenwich, ISS, t0, t1, 'lower')
    all_times = concatenate([upper_times.tt, lower_times.tt])
    all_times.sort()
    assert (times.tt == all_times).all()

    # Check that the found times produce the correct data
    f = partial(_satellite_alt, greenwich.positives[-1], ISS.positives[-1])
    assert is_extreme(f, times.tt, 2*ms)
    
    # Check that same result is found if plain Topos is used as observer
    times2 = culminations(greenwich.positives[-1], ISS, t0, t1, 'all')
    assert (times.tt == times2.tt).all()


if __name__ == '__main__':
    import pytest
    pytest.main(['-s'])