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
