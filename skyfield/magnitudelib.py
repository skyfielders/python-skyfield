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
    ap_mag = -0.613 + distance_mag_factor + ph_ang_factor
    return ap_mag

if __name__ == '__main__':
    print(mercury_magnitude(0.310295423552, 1.32182643625754, 1.1677))
    print(mercury_magnitude(0.413629222334, 0.9264480871861, 90.1662))
