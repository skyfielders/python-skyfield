"""Test that Skyfield implements NOVAS-compatible Earth deflection.

Run `design/measure_earth_deflection.py` to draw a quick graph.

"""
from numpy import arange, diff
from skyfield import relativity
from skyfield.api import load, wgs84
from skyfield.framelib import ecliptic_frame as ef

# The expected offsets for deflection in the following tests have NOT
# been confirmed against an external authority; they are simply what
# Skyfield 1.51 did when these tests were written.  But at least they
# stop us from accidentally changing Skyfield's behavior going forward.
#
# It's encouraging that Skyfield thinks that Earth deflection affects
# altitude but not azimuth, and bumps altitude up.

def stringify(angle1, angle2):
    return ' '.join('%.2f' % n for n in angle1.mas() - angle2.mas())

def test_sun_deflection():
    ts = load.timescale()
    t = ts.tt(2015, 3, 16, arange(0, 25, 6))
    eph = load('de421.bsp')
    earth = eph['earth']
    jupiter = eph['jupiter barycenter']
    lowell = earth + wgs84.latlon(35.2029, -111.6646)
    lat1, lon1, _ = lowell.at(t).observe(jupiter).apparent().frame_latlon(ef)

    old = relativity.rmasses['sun']
    relativity.rmasses['sun'] = 1e100  # almost no mass
    lat2, lon2, _ = lowell.at(t).observe(jupiter).apparent().frame_latlon(ef)
    relativity.rmasses['sun'] = old

    assert stringify(lat1, lat2) == '-0.02 -0.02 -0.02 -0.02 -0.02'
    assert stringify(lon1, lon2) == '1.27 1.28 1.29 1.30 1.31'

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

    # These values have not been checked against an external authority.
    # They are simply what Skyfield 1.51 did when tested.  But at least
    # this test prevents us from changing behavior in the future without
    # knowing it.  It's at least encouraging that Earth deflection
    # affects altitude not azimuth, and bumps altitude up.

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
