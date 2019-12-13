import numpy as np
from skyfield.api import PlanetaryConstants, T0, load
from skyfield.positionlib import ICRF

def test_frame_rotation():
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
    delta = r - spiceypy_matrix
    assert (delta < 1e-16).all()  # we agree to roughly float64 precision!

    # Second, a moment when the angle W is more than 2500 radians.

    tdb = T0
    spiceypy_matrix = [
        [ 0.7840447406961362,  0.5582359944893811,  0.2713787372716964],
        [-0.6203032939745002,  0.7203957219351799,  0.3102480093439375],
        [-0.0223084753202375, -0.4115854446818337,  0.9110981032001678],
    ]
    r = frame.rotation_at(ts.tdb_jd(tdb))
    delta = r - spiceypy_matrix
    assert (delta < 2e-13).all()  # 4 digits are lost in large W radians

def test_frame_rotation2():
    et_seconds = 259056665.0
    ts = load.timescale(builtin=True)
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
    print(delta)
    assert (delta < tolerance).all()

    # Example of from moon_080317.tf (hoping that the vector quoted below
    # is exactly the same one that is computed behind the scenes in
    # their example):

    vector = ICRF([247982.584262823, -261900.0687693838, -124558.3708890628])
    vector.t = t
    result = vector.frame_xyz(frame)
    print(result.au)
    assert max(abs(result.au - [379908.634, 33385.003, -12516.8859])) < 1e-1

    # TODO: Why is the agreement not nearly 1e-4?  It would be nice to
    # agree with the moon_080317.tf numbers to their level of precision,
    # which is only 9 digits.

    relative_frame = pc.build_frame_named('MOON_ME_DE421')
    result = vector.frame_xyz(relative_frame)
    print(result.au)
    print(result.au - [379892.825, 33510.118, -12661.5278])
    assert max(abs(result.au - [379892.825, 33510.118, -12661.5278])) < 1e-1

    # TODO:
    # a Moon-based topos object, tested against HORIZONS examples in repository
    # sensitive test of rotation based frame.

    # pc.read_text(load('pck00008.tpc'))
    # print(pc.assignments['BODY301_RADII'])

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
