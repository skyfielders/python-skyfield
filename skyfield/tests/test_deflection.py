"""Test that Skyfield implements NOVAS-compatible Earth deflection.

Run `design/measure_earth_deflection.py` to draw a quick graph.

"""
from numpy import arange, diff
from skyfield import relativity
from skyfield.api import load, wgs84
from skyfield.framelib import ecliptic_frame as ef
from skyfield.jpllib import _jpl_name_code_dict as jpl_codes

def trial(deflector_name, target_name, year, month, day):
    """Compare apparent() lat/lon with a given deflector turned on and off."""
    ts = load.timescale()
    t = ts.tt(year, month, day, arange(0, 25, 6))
    eph = load('de421.bsp')
    earth = eph['earth']
    target = eph[target_name]
    lowell = earth + wgs84.latlon(35.2029, -111.6646)
    astrometric = lowell.at(t).observe(target)
    lat1, lon1, _ = astrometric.apparent().frame_latlon(ef)

    code = jpl_codes[deflector_name.upper()]
    old = relativity.rmasses[code]
    relativity.rmasses[code] = 1e100  # almost no mass
    lat2, lon2, _ = astrometric.apparent().frame_latlon(ef)
    relativity.rmasses[code] = old

    return lat1.mas() - lat2.mas(), lon1.mas() - lon2.mas()

def stringify(difference_mas):
    return ' '.join('%.2f' % n for n in difference_mas)

# The expected offsets for deflection in the following tests have NOT
# been confirmed against an external authority; they are simply what
# Skyfield 1.51 did when these tests were written.  But at least they
# stop us from accidentally changing Skyfield's behavior going forward.

def test_sun_deflection():
    lat_mas, lon_mas = trial('sun', 'jupiter barycenter', 2025, 3, 16)
    assert stringify(lat_mas) == '-0.01 -0.01 -0.01 -0.01 -0.01'
    assert stringify(lon_mas) == '4.16 4.17 4.19 4.21 4.22'

def test_jupiter_deflection():
    # Very close conjunction of 0.1 degrees.
    lat_mas, lon_mas = trial('jupiter', 'saturn barycenter', 2020, 12, 21)
    assert stringify(lat_mas) == '0.20 0.25 0.31 0.33 0.31'
    assert stringify(lon_mas) == '0.16 0.14 0.09 0.01 -0.08'

def test_saturn_deflection():
    # Conjunction, fairly distant 0.83 degrees separation.
    lat_mas, lon_mas = trial('saturn', 'neptune barycenter', 2026, 2, 20)
    assert stringify(lat_mas) == '0.01 0.01 0.01 0.01 0.01'
    assert stringify(lon_mas) == '0.00 0.00 0.00 -0.00 -0.00'

def test_moon_deflection():  # TODO: better dates
    # Moon, the next deflector, should be turned off by default.
    lat_mas, lon_mas = trial('moon', 'neptune barycenter', 2024, 12, 9.364)
    assert stringify(lat_mas) == '0.00 0.00 0.00 0.00 0.00'
    assert stringify(lon_mas) == '0.00 0.00 0.00 0.00 0.00'

# Not tested against an external authority; but it's encouraging that
# Earth deflection affects altitude not azimuth, and bumps altitude up
# not down.

def test_earth_deflection():
    ts = load.timescale()
    t = ts.tt(2015, 3, 16, arange(0, 25, 2))
    eph = load('de421.bsp')
    earth = eph['earth']
    jupiter = eph['jupiter barycenter']
    lowell = earth + wgs84.latlon(35.2029, -111.6646)
    alt1, az1, _ = lowell.at(t).observe(jupiter).apparent().altaz()

    old = relativity.rmasses['earth']
    relativity.rmasses['earth'] = 1e100  # almost no mass
    alt2, az2, _ = lowell.at(t).observe(jupiter).apparent().altaz()
    relativity.rmasses['earth'] = old

    dalt = ' '.join('%.2f' % n for n in alt1.mas() - alt2.mas())
    daz = ' '.join('%.2f' % n for n in az1.mas() - az2.mas())

    s = '0.19 0.11 0.06 0.06 0.12 0.19 0.30 0.00 0.00 0.00 0.00 0.29 0.19'
    assert dalt == s

    s = '0.00 0.00 0.00 -0.00 -0.00 -0.00 -0.00 0.00 0.00 0.00 0.00 0.00 0.00'
    assert daz == s

def test_earth_deflection_cutoff():
    # The NOVAS library includes the Earth's gravitational deflection of
    # light for both topocentric observers and observers in Earth orbit,
    # but shuts this effect off once the object is behind the Earth 20%
    # of the way from its limb towards its center.  This test determines
    # whether Skyfield puts the resulting discontinuity in the same
    # place as the NOVAS library does.

    ts = load.timescale(delta_t=0.0)
    t = ts.tt(2016, 7, 2, arange(10.5628, 10.5639, 0.0002))
    planets = load('de405.bsp')
    earth = planets['earth']
    mars = planets['mars']
    lowell = earth + wgs84.latlon(35.2029, -111.6646)
    ra, dec, distance = lowell.at(t).observe(mars).apparent().radec()
    h = ra.hours
    hprime = diff(h)
    assert hprime[0] > 1.8e-8
    assert hprime[1] > 1.8e-8
    assert hprime[2] < 1.3e-8   # moment when nadir angle crosses 0.8
    assert hprime[3] > 1.8e-8
    assert hprime[4] > 1.8e-8
