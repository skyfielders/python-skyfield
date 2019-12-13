import numpy as np
from skyfield.api import PlanetaryConstants, T0, load
from skyfield.positionlib import ICRF

def test_frame_rotation_matrices():
    # Ask a frame defined directly by a binary segment for rotation
    # matrices, to test how we build and compose the rotation matrices.

    ts = load.timescale(builtin=True)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))
    frame = pc.build_frame_named('MOON_PA_DE421')
    assert frame._matrix is None  # pure segment, with no rotation applied

    # First, a moment when the angle W is nearly zero radians, so all of
    # its precision is available to the trigonometry.

    tdb = T0 - 11150
    spiceypy_matrix = [
        [ 0.9994150897380264,  0.0323102706039267,  0.0112037858527199],
        [-0.0341574268117634,  0.9272642685944782,  0.3728461430423364],
        [ 0.0016578894811168, -0.3730107540024127,  0.9278255379116378],
    ]
    r = frame.rotation_at(ts.tdb_jd(tdb))
    delta = abs(r - spiceypy_matrix)
    assert (delta < 6e-17).all()  # nearly float64 precision

    # Second, a moment when the angle W is more than 2500 radians.

    tdb = T0
    spiceypy_matrix = [
        [ 0.7840447406961362,  0.5582359944893811,  0.2713787372716964],
        [-0.6203032939745002,  0.7203957219351799,  0.3102480093439375],
        [-0.0223084753202375, -0.4115854446818337,  0.9110981032001678],
    ]
    r = frame.rotation_at(ts.tdb_jd(tdb))
    delta = abs(r - spiceypy_matrix)
    assert (delta < 2e-13).all()  # a few digits are lost in large W radians

def test_rotating_vector_into_frame():
    et_seconds = 259056665.1855896
    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(T0 + et_seconds / 3600. / 24.0)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    # Example from "moon_080317.tf" (the raw vector is taken from a
    # tweaked version of the script in the file that uses "J2000"):

    vector = ICRF([2.4798273371071659e+05,
                   -2.6189996683651494e+05,
                   -1.2455830876097400e+05])
    vector.t = t

    frame = pc.build_frame_named('MOON_PA_DE421')
    result = vector.frame_xyz(frame)
    assert max(abs(result.au - [379908.634, 33385.003, -12516.8859])) < 1e-3

    relative_frame = pc.build_frame_named('MOON_ME_DE421')
    result = vector.frame_xyz(relative_frame)
    assert max(abs(result.au - [379892.825, 33510.118, -12661.5278])) < 1e-2

    # TODO:
    # a Moon-based topos object, tested against HORIZONS examples in repository
    # sensitive test of rotation based frame.

def test_position_of_latitude_longitude_on_moon():
    ts = load.timescale(builtin=True)
    t = ts.utc(2019, 12, 10)  # TODO: switch to TDB date?
    print(t.tdb)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    from ..constants import AU_KM, tau
    from ..functions import rot_z, rot_y

    print(360 - 46.8)
    lat = 26.3 * tau / 360.0  # north
    lon = -46.8 * tau / 360.0  # east (negative: west)

    frame = pc.build_frame_named('MOON_ME_DE421')
    position = np.array((1737.4, 0, 0))
    position = rot_y(-lat).dot(position)
    position = rot_z(lon).dot(position)  # sign?

    from ..functions import length_of
    print('length km:', length_of(np.array((1.043588965592271E-05, 1.534286674325492E-06, -4.859893068052085E-06))) * AU_KM)

    print(position)

    position /=  AU_KM

    R = frame.rotation_at(t)
    #position = R.dot(position)
    position = R.T.dot(position)
    position *= -1.0
    print(position)
    print(position-[1.043588965592271E-05, 3.340834944508400E-06, -3.848560523814720E-06])

    want = [1.043588965592271E-05, 3.340834944508400E-06,
            -3.848560523814720E-06]
    assert max(abs(position - want)) < 2e-9

    # from ..planetarylib import PlanetTopos
    # topos = PlanetTopos(301, frame, position)
    # geometric = topos.at(t)

    # print(geometric.position.au)

def test_observing_from_location_on_moon():
    ts = load.timescale(builtin=True)
    t = ts.utc(2019, 12, 13)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    from skyfield.constants import AU_KM
    frame = pc.build_frame_named('MOON_ME_DE421')
    position = np.array((1737.1, 0, 0)) / AU_KM

    from ..planetarylib import PlanetTopos
    topos = PlanetTopos(301, frame, position)

    print(topos.at(t))

def test_frame_alias():
    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    f1 = pc.build_frame_named('MOON_PA_DE421')
    f2 = pc.build_frame_named('MOON_PA')

    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(T0)
    assert (f1.rotation_at(t) == f2.rotation_at(t)).all()

def test_unloaded_bpc():
    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    try:
        pc.build_frame_named('MOON_PA_DE421')
    except LookupError as e:
        assert str(e) == (
            'you have not yet loaded a binary PCK file'
            ' that has a segment for frame 31006'
        )
