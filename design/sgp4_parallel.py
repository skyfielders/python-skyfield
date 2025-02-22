#!/usr/bin/env python3
"""
The SGP4 library has the ability to process an array of satellites across an
array of times. What could support for this look like in Skyfield, combining
optimized parallel calculation (in the C++ version of SGP4) with coordinate
conversions from TEME?
"""
import numpy as np
from sgp4.api import SatrecArray, SGP4_ERRORS, accelerated
from skyfield.api import load, utc, EarthSatellite, Topos
from skyfield.constants import ANGVEL, AU_KM, DAY_S, AU_M, tau
from skyfield.sgp4lib import TEME_to_ITRF, _zero_zero_minus_one, theta_GMST1982, _cross
from skyfield.positionlib import build_position
from skyfield.vectorlib import ObserverData
from skyfield.functions import mxv, rot_z

# Precompute values to be reused across all satellites for the same times for conversion from TEME
def precompute_for_TEME(jd_ut1):
    theta, theta_dot = theta_GMST1982(jd_ut1)
    angular_velocity = np.multiply.outer(_zero_zero_minus_one, theta_dot)
    R = rot_z(-theta)
    return angular_velocity, R

# Based on an earlier version of TEME_to_ITRF() and does not include xp, yp, or fraction_ut1
def TEME_to_ITRF_fast(jd_ut1, rTEME, vTEME, angular_velocity, R):
    if len(rTEME.shape) == 1:
        rPEF = (R).dot(rTEME)
        vPEF = (R).dot(vTEME) + _cross(angular_velocity, rPEF)
    else:
        rPEF = mxv(R, rTEME)
        vPEF = mxv(R, vTEME) + _cross(angular_velocity, rPEF)
    return rPEF, vPEF

def precompute_for_ITRF(t):
    spin = rot_z(t.gast / 24.0 * tau)
    return spin

# Based on an earlier version of ITRF_to_GCRS2() and does not include _high_accuracy for now
def ITRF_to_GCRS2_fast(t, rITRF, vITRF, spin):
    position = mxv(spin, rITRF)
    velocity = mxv(spin, vITRF)
    velocity[0] += DAY_S * ANGVEL * - position[1]
    velocity[1] += DAY_S * ANGVEL * position[0]
    position = mxv(t.MT, position)
    velocity = mxv(t.MT, velocity)
    return position, velocity

# Run SGP4 in parallel across n satellites x m times
def positions_for(satellites, earthLocation, times):
    sat_array = SatrecArray( [s.model for s in satellites] )
    jd = times._utc_float()
    e, r, v = sat_array.sgp4(jd, np.zeros_like(jd))
    # Cached computations
    _, loc_v_GCRS, loc_p_GCRS, _ = earthLocation._at(times)
    loc_altaz_rotation = earthLocation._altaz_rotation(times)
    no_errors = [None] * len(e[0])
    angular_velocity, R = precompute_for_TEME(times.ut1)
    spin = precompute_for_ITRF(times)
    # Unpack the TEME coordinates, convert to GCRS and subtract earthLocation
    for index in range(r.shape[0]):
        errors = e[index]
        messages = no_errors if not np.any(errors) else \
                    [SGP4_ERRORS[error] if error else None for error in errors ]
        # Adapted from _position_and_velocity_TEME_km(), ITRF_position_velocity_error()
        rTEME = np.divide(r[index].T, AU_KM)
        vTEME = np.divide(v[index].T, AU_KM / DAY_S)
        rITRF, vITRF = TEME_to_ITRF_fast(times.ut1, rTEME, vTEME, angular_velocity, R)
        sat_p_GCRS, sat_v_GCRS = ITRF_to_GCRS2_fast(times, rITRF, vITRF, spin)
        #sat_p_GCRS, sat_v_GCRS = ITRF_to_GCRS2(times, rITRF, vITRF)
        # Mirror VectorSum _at() for (EarthSatellite - Topos)
        sat_p_GCRS -= loc_p_GCRS
        sat_v_GCRS -= loc_v_GCRS
        # Mirror VectorFunction at()
        observer_data = ObserverData()
        observer_data.gcrs_position = loc_p_GCRS
        observer_data.altaz_rotation = loc_altaz_rotation
        observer_data.elevation_m = earthLocation.elevation.m
        target = -100000 - satellites[index].model.satnum
        position = build_position(sat_p_GCRS, sat_v_GCRS, times, earthLocation, target, observer_data)
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
    satellites = [iss]
    city = Topos('44.0247 N', '88.5426 W')
    for geo1 in positions_for(satellites, city, times):
        print(geo1.xyz.km)
    #
    print("Compare to EarthSatellite at():")
    difference = iss - city
    geo2 = difference.at(times)
    #geo2 = iss.at(times)
    print(geo2.xyz.km)
    #
    assert np.allclose(geo1.xyz.km, geo2.xyz.km)
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
    city = Topos('44.0247 N', '88.5426 W')
    trange = range(0, 300, 10)
    times = ts.utc(2020, 7, 5, 12, 0, trange)
    satellites = [kestrel, iss]
    positions = list(positions_for(satellites, city, times))
    kestrel1, iss1 = positions
    iss2 = (iss - city).at(times)
    kestrel2 = (kestrel - city).at(times)
    assert np.allclose(iss1.xyz.km, iss2.xyz.km)
    assert np.allclose(kestrel1.xyz.km, kestrel2.xyz.km)
    assert np.allclose(iss1.velocity.au_per_d, iss2.velocity.au_per_d)
    assert np.allclose(kestrel1.velocity.au_per_d, kestrel2.velocity.au_per_d)

def test_many_sats():
    print("test many sats")
    ts = load.timescale(builtin=True)
    active_sats = load.tle_file("https://celestrak.org/NORAD/elements/active.txt", reload=False)
    print(" active satellites:", len(active_sats), " x 300 time values")
    trange = range(0, 300, 1)
    times = ts.utc(2020, 8, 25, 12, 0, trange)
    city = Topos('44.0247 N', '88.5426 W')
    positions = list(positions_for(active_sats, city, times))
    assert len(active_sats) == len(positions)

if __name__ == "__main__":
    main()
