"""Hand-crafted tests against specific NOVAS behaviors.

The tests in the neighboring `tests_against_novas.py` are automatically
generated and in general demonstrate close agreement with NOVAS.  But
the hand-crafted tests here are aimed at specific edge conditions that
we want to make sure we get correct.

"""
from numpy import arange, diff
from skyfield.api import Topos, load

def test_earth_deflection():
    # The NOVAS library includes the Earth's gravitational deflection of
    # light for both topocentric observers and observers in Earth orbit,
    # but shuts this effect off once the object is behind the Earth 20%
    # of the way from its limb towards its center.  This test determines
    # whether Skyfield puts the resulting discontinuity in the same
    # place as the NOVAS library does.
    #
    # For more details see:
    # https://github.com/skyfielders/astronomy-notebooks/blob/master/Skyfield-Notes/Fixing-earth-deflection.ipynb

    t = load.timescale(delta_t=0.0)
    t = t.tt(2016, 7, 2, arange(10.5628, 10.5639, 0.0002))
    planets = load('de405.bsp')
    earth = planets['earth']
    mars = planets['mars']
    lowell = earth + Topos(latitude_degrees=35.2029, longitude_degrees=-111.6646)
    ra, dec, distance = lowell.at(t).observe(mars).apparent().radec()
    h = ra.hours
    hprime = diff(h)
    assert hprime[0] > 1.8e-8
    assert hprime[1] > 1.8e-8
    assert hprime[2] < 1.3e-8   # moment when nadir angle crosses 0.8
    assert hprime[3] > 1.8e-8
    assert hprime[4] > 1.8e-8
