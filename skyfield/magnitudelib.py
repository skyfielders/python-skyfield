"""Routines for computing magnitudes.

Planetary routines adapted from:

https://arxiv.org/pdf/1808.01973.pdf

"""
from numpy import log10

def mercury_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    ph_ang_factor = (
        6.3280e-02 * ph_ang
        - 1.6336e-03 * ph_ang**2
        + 3.3644e-05 * ph_ang**3
        - 3.4265e-07 * ph_ang**4
        + 1.6893e-09 * ph_ang**5
        - 3.0334e-12 * ph_ang**6
    )
    return -0.613 + distance_mag_factor + ph_ang_factor

def venus_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    if ph_ang < 163.7:
        ph_ang_factor = (
            -1.044E-03 * ph_ang
            + 3.687E-04 * ph_ang**2
            - 2.814E-06 * ph_ang**3
            + 8.938E-09 * ph_ang**4
        )
    else:
        ph_ang_factor = (
            236.05828 + 4.384
            - 2.81914E+00 * ph_ang
            + 8.39034E-03 * ph_ang**2
        )
    return -4.384 + distance_mag_factor + ph_ang_factor

def earth_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10 (r * delta)
    ph_ang_factor = -1.060e-03 * ph_ang + 2.054e-04 * ph_ang**2
    return -3.99 + distance_mag_factor + ph_ang_factor

def jupiter_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    geocentric_phase_angle_limit = 12.0

    if ph_ang <= geocentric_phase_angle_limit:
        ph_ang_factor = -3.7E-04 * ph_ang + 6.16E-04 * ph_ang**2
    else:
        ph_ang_factor = -2.5 * log10(
            1.0 - 1.507 * (ph_ang / 180.)
            - 0.363 * (ph_ang / 180.)**2
            - 0.062 * (ph_ang / 180.)**3
            + 2.809 * (ph_ang / 180.)**4
            - 1.876 * (ph_ang / 180.)**5
        )

    if ph_ang <= geocentric_phase_angle_limit:
        ap_mag = -9.395 + distance_mag_factor + ph_ang_factor
    else:
        ap_mag = -9.428 + distance_mag_factor + ph_ang_factor

    return ap_mag

def uranus_magnitude(r, delta, ph_ang,
                     sun_sub_lat_planetog, earth_sub_lat_planetog):
    distance_mag_factor = 5.0 * log10 (r * delta)
    sub_lat_planetog = (abs(sun_sub_lat_planetog)
                        + abs(earth_sub_lat_planetog)) / 2.0
    sub_lat_factor = -0.00084 * sub_lat_planetog
    geocentric_phase_angle_limit = 3.1
    ap_mag = -7.110 + distance_mag_factor + sub_lat_factor
    if ph_ang > geocentric_phase_angle_limit:
        ap_mag += 6.587e-3 * ph_ang + 1.045e-4 * ph_ang**2
    return ap_mag
