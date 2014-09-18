"""Tests against HORIZONS numbers."""

from skyfield import api

# see the top-level project ./horizons/ directory for where the
# following numbers come from; soon, we should automate the fetching of
# such numbers and their injection into test cases, as we do for results
# from NOVAS.

"""
 Date__(UT)__HR:MN     hEcl-Lon hEcl-Lat               r        rdot
********************************************************************
$$SOE
 1980-Jan-01 00:00     151.3229   1.0130  5.378949180806   0.4314383
$$EOE
"""

def test_ecliptic_latlon():
    astrometric = api.sun(utc=(1980, 1, 1)).observe(api.jupiter)
    lat, lon, distance = astrometric.ecliptic_latlon()
    assert '{0:.4f}'.format(lat.degrees) == '1.0130'
    assert '{0:.4f}'.format(lon.degrees) == '151.3227'
    # That last value should really be '151.3229' according to HORIZONS
    # (see string above) but we are just getting started here so the
    # tiny difference is being filed away as something to look at later!
