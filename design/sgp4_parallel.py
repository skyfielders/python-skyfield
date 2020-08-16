#!/usr/bin/env python3
"""
The SGP4 library has the ability to process an array of satellites across an
array of times. What could support for this look like in Skyfield, combining
optimized parallel calculation (in the C++ version of SGP4) with coordinate
conversions from TEME?
"""
import numpy as np
from sgp4.api import SatrecArray, SGP4_ERRORS, accelerated
from skyfield.api import load, utc, EarthSatellite
from skyfield.constants import AU_KM, DAY_S, AU_M
from skyfield.sgp4lib import TEME_to_ITRF
from skyfield.positionlib import ITRF_to_GCRS2, build_position
from skyfield.vectorlib import ObserverData

# SGP4 in parallel across n satellites x m times
def positions_for(satellites, times):
    sat_array = SatrecArray( [s.model for s in satellites] )
    jd = times._utc_float()
    e, r, v = sat_array.sgp4(jd, np.zeros_like(jd))
    assert len(satellites) == np.shape(r)[0]
    for index in range(r.shape[0]):
        errors = e[index]
        messages = [SGP4_ERRORS[error] if error else None for error in errors ]
        # Adapted from _position_and_velocity_TEME_km(), ITRF_position_velocity_error()
        rTEME = np.divide(r[index].T, AU_KM)
        vTEME = np.divide(v[index].T, AU_KM / DAY_S)
        rITRF, vITRF = TEME_to_ITRF(times.ut1, rTEME, vTEME)
        rGCRS, vGCRS = ITRF_to_GCRS2(times, rITRF, vITRF)
        # Adapted from VectorFunction at()
        observer_data = ObserverData()
        observer_data.gcrs_position = rGCRS
        target = -100000 - satellites[index].model.satnum
        position = build_position(rGCRS, vGCRS, times, 399, target, observer_data)
        position.message = messages
        yield position

# ------------------------------------------------------------------------- #

ISS_TLE = """1 25544U 98067A   20187.34541214  .00002866  00000-0  59793-4 0  9991
2 25544  51.6454 258.5877 0002990  76.8215 328.5773 15.49181247234835
"""

KESTREL_TLE = """1 42982U 98067NE  20187.19650118  .00008987  00000-0  76592-4 0  9999
2 42982  51.6323 190.7341 0002835 141.7806 218.3393 15.70424148153865
"""

def main():
    if not accelerated:
        print("SGP4 library is not accelerated")
        return
    test_iss()
    test_two_sats()
    test_many_sats()

def test_iss():
    ts = load.timescale(builtin=True)
    L1, L2 = ISS_TLE.splitlines()
    iss = EarthSatellite(L1, L2)

    trange = range(0, 60, 10)
    times = ts.utc(2020, 7, 5, 12, 0, trange)
    print("positions_for():")
    satellites = [iss, iss]
    for geo1 in positions_for(satellites, times):
        print(geo1.position.km)
        #print(geo1.message)
    #
    print("Compare to EarthSatellite at():")
    geo2 = iss.at(times)
    print(geo2.position.km)
    #
    # Tiny velocity difference compared to `at()` triggered by moving `AU_KM / DAY_S` into np.divide()
    #np.set_printoptions(precision=16)
    #print(geo1.velocity.km_per_s - geo2.velocity.km_per_s)
    #
    assert np.isclose(geo1.position.km, geo2.position.km).all()
    one_m_per_hour = 1.0 * 24.0 / AU_M
    assert abs(geo1.velocity.au_per_d - geo2.velocity.au_per_d).max() < one_m_per_hour

def test_two_sats():
    print("test two sats")
    ts = load.timescale(builtin=True)
    L1, L2 = ISS_TLE.splitlines()
    iss = EarthSatellite(L1, L2)
    L1, L2 = KESTREL_TLE.splitlines()
    kestrel = EarthSatellite(L1, L2)
    #
    trange = range(0, 60, 10)
    times = ts.utc(2020, 7, 5, 12, 0, trange)
    satellites = [kestrel, iss]
    positions = list(positions_for(satellites, times))
    kestrel1, iss1  = positions
    iss2 = iss.at(times)
    kestrel2 = kestrel.at(times)
    assert np.isclose(iss1.position.km, iss2.position.km).all()
    assert np.isclose(kestrel1.position.km, kestrel2.position.km).all()

# With optimized SGP4, peak memory use 142MB and runs in ~1s on my laptop
def test_many_sats():
    print("test many sats")
    ts = load.timescale(builtin=True)
    active_sats = load.tle_file("https://celestrak.com/NORAD/elements/active.txt", reload=False)
    print("  #active satellites:", len(active_sats), " x 300 time values")
    trange = range(0, 300, 1)
    times = ts.utc(2020, 7, 5, 12, 0, trange)
    locations = list(positions_for(active_sats, times))
    assert len(active_sats) == len(locations)

if __name__ == "__main__":
    main()
