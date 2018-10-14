from skyfield.api import Loader, Topos, Angle, EarthSatellite, Star
from skyfield.timelib import Time
from almanac2 import (meridian_transits, culminations, twilights, 
                       risings_settings, seasons, moon_phases)
from almanac2 import (_ecliptic_lon_diff, _moon_ul_alt, _lha, _alt,
                       _satellite_alt, _ecliptic_lon)
from numpy import ndarray
from functools import partial

load = Loader(r'C:\Users\Josh\Scripts\Skyfield_Data')

ts = load.timescale()
ephem = load('de430.bsp')
earth = ephem['earth']
sun = ephem['sun']
moon = ephem['moon']
mars = ephem['mars barycenter']

greenwich = Topos('51.5 N', '0 W')

iss_tle = """\
1 25544U 98067A   18161.85073725  .00003008  00000-0  52601-4 0  9993
2 25544  51.6418  50.3007 0003338 171.6979 280.7366 15.54163173117534
"""
ISS = EarthSatellite(*iss_tle.splitlines())

sirius = Star(ra_hours=(6, 45, 8.91728), dec_degrees=(-16, 42, 58.0171))

minute = 1/24/60
ms = minute/60/1000


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
    times, lons = seasons(earth, t0, t1)
    
    assert ((lons.degrees==0) + (lons.degrees==90) + (lons.degrees==180) + (lons.degrees==270)).all()
    
    assert isinstance(times, Time)
    assert isinstance(lons, Angle)
    assert len(times) == len(lons.radians) == 104
    
    compare(times[0].tt, ts.utc(2000, 3, 20, 7, 35).tt, minute/2)
    compare(times[1].tt, ts.utc(2000, 6, 21, 1, 48).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_ecliptic_lon, earth, sun)
    assert is_root(f, times.tt, lons.degrees, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_phases.pl?year=2017&month=1&day=1&nump=50&format=p
def test_moon_phases():
    t0 = ts.utc(2017)
    t1 = ts.utc(2018)
    times, lon_diffs = moon_phases(moon, t0, t1)
    
    assert ((lon_diffs.degrees==0) + (lon_diffs.degrees==90) + (lon_diffs.degrees==180) + (lon_diffs.degrees==270)).all()

    assert isinstance(times, Time)
    assert isinstance(lon_diffs, Angle)
    assert len(times) == len(lon_diffs.radians) == 49
    
    compare(times[0].tt, ts.utc(2017, 1, 5, 19, 47).tt, minute/2)
    
    # Check that the found times produce the correct data    
    f = partial(_ecliptic_lon_diff, earth, moon, sun)
    assert is_root(f, times.tt, lon_diffs.degrees, ms/2)


# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=0&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_sun_risings_settings():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, kinds = risings_settings(greenwich, sun, t0, t1)
    
    assert ((kinds=='rise') + (kinds=='set')).all()
    
    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 8, 6).tt, minute/2)
    
    # Check that the found times produce the correct data    
    f = partial(_alt, earth+greenwich, sun)
    assert is_root(f, times.tt, -50/60, ms/2)


# Data Source:   
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=1&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_moon_risings_settings():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, kinds = risings_settings(greenwich, moon, t0, t1)
    
    assert ((kinds=='rise') + (kinds=='set')).all()
    
    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 60
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 9, 46).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_moon_ul_alt, earth+greenwich, moon)
    assert is_root(f, times.tt, -34/60, ms/2)
    
    
def test_ISS_risings_settings():
    t0 = ts.utc(2017, 6, 1)
    t1 = ts.utc(2017, 6, 2)
    times, kinds = risings_settings(greenwich, ISS, t0, t1)
    
    assert ((kinds=='rise') + (kinds=='set')).all()
    
    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 13
    
    # Check that the found times produce the correct data
    f = partial(_satellite_alt, greenwich, ISS)
    assert is_root(f, times.tt, -34/60, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=2&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_civil_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, am_pm = twilights(greenwich, sun, t0, t1, 'civil')

    assert ((am_pm=='am') + (am_pm=='pm')).all()

    assert isinstance(times, Time)
    assert isinstance(am_pm, ndarray)
    assert len(times) == len(am_pm) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 7, 26).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_alt, earth+greenwich, sun)
    assert is_root(f, times.tt, -6, ms/2)


# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=3&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_nautical_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, am_pm = twilights(greenwich, sun, t0, t1, 'nautical')
    
    assert ((am_pm=='am') + (am_pm=='pm')).all()

    assert isinstance(times, Time)
    assert isinstance(am_pm, ndarray)
    assert len(times) == len(am_pm) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 6, 43).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_alt, earth+greenwich, sun)
    assert is_root(f, times.tt, -12, ms/2)
    
    
# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=4&place=Greenwich&lon_sign=-1&lon_deg=0&lon_min=0&lat_sign=1&lat_deg=51&lat_min=30&tz=&tz_sign=-1
def test_astronomical_twilights():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, am_pm = twilights(greenwich, sun, t0, t1, 'astronomical')
    
    assert ((am_pm=='am') + (am_pm=='pm')).all()

    assert isinstance(times, Time)
    assert isinstance(am_pm, ndarray)
    assert len(times) == len(am_pm) == 62
    
    compare(times[0].tt, ts.utc(2017, 1, 1, 6, 2).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_alt, earth+greenwich, sun)
    assert is_root(f, times.tt, -18, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_mrst2.pl?form=2&ID=AA&year=2017&month=1&day=1&reps=31&body=10&place=Greenwich&lon_sign=-1&lon_deg=&lon_min=&lon_sec=&lat_sign=1&lat_deg=51&lat_min=30&lat_sec=&height=&tz=&tz_sign=-1
def test_sun_transits():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, lhas = meridian_transits(greenwich, sun, t0, t1)
    
    assert ((lhas.hours==0) + (lhas.hours==12)).all()

    assert isinstance(times, Time)
    assert isinstance(lhas, Angle)
    assert len(times) == len(lhas.radians) == 62
        
    compare(times[1].tt, ts.utc(2017, 1, 1, 12, 4).tt, minute/2)

    # Check that the found times produce the correct data
    f = partial(_lha, earth+greenwich, sun)
    assert is_root(f, times.tt, lhas._degrees, ms/2)
    

# Data Source:
# http://aa.usno.navy.mil/cgi-bin/aa_mrst2.pl?form=2&ID=AA&year=2017&month=1&day=1&reps=31&body=4&place=Greenwich&lon_sign=-1&lon_deg=&lon_min=&lon_sec=&lat_sign=1&lat_deg=51&lat_min=30&lat_sec=&height=&tz=&tz_sign=-1
def test_mars_transits():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, lhas = meridian_transits(greenwich, mars, t0, t1)
    
    assert ((lhas.hours==0) + (lhas.hours==12)).all()
    
    assert isinstance(times, Time)
    assert isinstance(lhas, Angle)
    assert len(times) == len(lhas.radians) == 62
    
    compare(times[1].tt, ts.utc(2017, 1, 1, 16, 2).tt, minute/2)
    
    # Check that the found times produce the correct data
    f = partial(_lha, earth+greenwich, mars)
    assert is_root(f, times.tt, lhas._degrees, ms/2)
    

def test_sun_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, kinds = culminations(greenwich, sun, t0, t1)
    
    assert ((kinds=='upper') + (kinds=='lower')).all()
    
    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 62
    
    # Check that the found times produce the correct data
    f = partial(_alt, earth+greenwich, sun)        
    assert is_extreme(f, times.tt, 2*ms)
    
    
def test_moon_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 32)
    times, kinds = culminations(greenwich, moon, t0, t1)
    
    assert ((kinds=='upper') + (kinds=='lower')).all()

    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 60

    # Check that the found times produce the correct data
    f = partial(_alt, earth+greenwich, moon)    
    assert is_extreme(f, times.tt, 6*ms)
    
    
def test_ISS_culminations():
    t0 = ts.utc(2017, 1, 1)
    t1 = ts.utc(2017, 1, 2)
    times, kinds = culminations(greenwich, ISS, t0, t1)
    
    assert ((kinds=='upper') + (kinds=='lower')).all()

    assert isinstance(times, Time)
    assert isinstance(kinds, ndarray)
    assert len(times) == len(kinds) == 31

    # Check that the found times produce the correct data
    f = partial(_satellite_alt, greenwich, ISS)
    assert is_extreme(f, times.tt, 2*ms)


def test_small_intervals():
    t0 = ts.utc(2018)
    t1 = ts.utc(2018, 1, 1, 1, 1)
    
    times, kinds = meridian_transits(greenwich, moon, t0, t1)
    assert len(times) == len(kinds.radians) == 0

    times, kinds = culminations(greenwich, moon, t0, t1)
    assert len(times) == len(kinds) == 0
    
    times, kinds = risings_settings(greenwich, moon, t0, t1)
    assert len(times) == len(kinds) == 0
    
    times, kinds = twilights(greenwich, sun, t0, t1)
    assert len(times) == len(kinds) == 0
    
    times, kinds = seasons(earth, t0, t1)
    assert len(times) == len(kinds.radians) == 0
    
    times, kinds = moon_phases(moon, t0, t1)
    assert len(times) == len(kinds.radians) == 0


if __name__ == '__main__':
    import pytest
    pytest.main(['-s'])