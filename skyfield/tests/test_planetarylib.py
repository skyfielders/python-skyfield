import numpy as np
from skyfield.api import PlanetaryConstants, T0, load

def test_rotation():
    et_seconds = 259056665.0
    ts = load.timescale()
    t = ts.tdb_jd(T0 + et_seconds / 3600. / 24.0)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.frame_from_name('MOON_PA_DE421')
    r = frame.rotation_at(t)

    delta = r - [
        [0.581793642428, -0.74958441165, -0.315657040855],
        [0.813188534322, 0.543494230073, 0.208178840241],
        [0.015510166907, -0.377805812142, 0.925754936814],
    ]
    assert (delta < 1e-10).all()  # TODO: why not 1e-15?
