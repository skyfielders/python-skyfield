import numpy as np
from skyfield.api import PlanetaryConstants, T0, load

def test_rotation():
    et_seconds = 259056665.0
    ts = load.timescale()
    t = ts.tdb_jd(T0 + et_seconds / 3600. / 24.0)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_PA_DE421')
    r = frame.rotation_at(t)

    # TODO: Why is this not better?  The rotation angle W is more than
    # 3,000 by the date we are testing, which throws away around 1e3
    # precision that a radian value might otherwise be expected to have.
    # But where does the other 1e3 go between the 1e16 precision of a
    # float64 and the actual precision of our matrix?
    tolerance = 1e-10

    # To generate the following matrix:
    #
    # import spiceypy as spice
    # spice.furnsh('moon_080317.tf')
    # spice.furnsh('moon_pa_de421_1900-2050.bpc')
    # spice.pxform('J2000', 'MOON_PA', et_seconds)
    spiceypy_matrix = [
        [0.5817936424278365, -0.7495844116502467, -0.3156570408552677],
        [0.8131885343221845, 0.5434942300733778, 0.2081788402405092],
        [0.0155101669071617, -0.3778058121415016, 0.9257549368135243],
    ]
    delta = r - spiceypy_matrix
    assert (delta < tolerance).all()
