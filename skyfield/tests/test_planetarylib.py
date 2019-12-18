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

    # To produce the following spiceypy matrices:
    #
    import numpy as np
    import spiceypy as spice
    from skyfield.constants import DAY_S
    spice.furnsh('moon_080317.tf')
    spice.furnsh('moon_pa_de421_1900-2050.bpc')
    g = np.vectorize(repr)
    print(g(spice.pxform('J2000', 'MOON_PA', (tdb - T0) * DAY_S)))
    print(g(spice.sxform('J2000', 'MOON_PA', (tdb - T0) * DAY_S)[3:6,:3]))

    spiceypy_rotation = [
        [0.9994150897380264, 0.032310270603926675, 0.011203785852719871],
        [-0.034157426811763446, 0.9272642685944782, 0.37284614304233643],
        [0.0016578894811167893, -0.3730107540024127, 0.9278255379116378],
    ]
    # spiceypy_rate = [
    #     [ -9.0916992653055378e-08,   2.4682369763316187e-06,
    #        9.9202268744491445e-07],
    #     [ -2.6601275013494018e-06,  -8.6208910075886084e-08,
    #       -2.9300741588246011e-08],
    #     [  4.2459631446005060e-10,  -5.0678787655029271e-10,
    #       -2.0450101227202991e-10],
    # ]
    spiceypy_rate = [
        [-9.091699265305538e-08, 2.4682369763316187e-06, 9.920226874449144e-07],
        [-2.660127501349402e-06, -8.620891007588608e-08, -2.930074158824601e-08],
        [4.245963144600506e-10, -5.067878765502927e-10, -2.045010122720299e-10],
    ]

    R = frame.rotation_at(ts.tdb_jd(tdb))
    assert (R == spiceypy_rotation).all()  # Boom.

    R2, Rv = frame.rotation_and_rate_at(ts.tdb_jd(tdb))
    assert (R == R2).all()  # (Less surprising: Python agrees with itself.)

    print('mine:')
    print(Rv)
    print('actual:')
    print(np.array(spiceypy_rate))
    print('difference:')
    print(Rv - spiceypy_rate)

    assert (Rv == spiceypy_rate).all()  # Boom.
    # print('candidate2:',
    #       repr(0 + -2.6616852779228065e-06 * 0.032310270603926675
    #            + -4.0568232474023717e-10 * 0.011203785852719871))
    print('candidate1:',
          repr(0 +  2.6616852779228065e-06 * -0.034157426811763446
               + 4.0568232474023717e-10 * 0.0016578894811167893))
    print('candidateC:',
          repr(0 +  2.6616852779228065e-06 *   -3.4157426811763446e-02
              + -4.0568232474023717e-10 *    1.6578894811167893e-03))

    # Second, a moment when the angle W is more than 2500 radians.

    tdb = T0
    spiceypy_rotation = [
        [ 0.7840447406961362,  0.5582359944893811,  0.2713787372716964],
        [-0.6203032939745002,  0.7203957219351799,  0.3102480093439375],
        [-0.0223084753202375, -0.4115854446818337,  0.9110981032001678],
    ]
    r = frame.rotation_at(ts.tdb_jd(tdb))
    delta = abs(r - spiceypy_rotation)
    assert (delta < 2e-13).all()  # a few digits are lost in large W radians?

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
    t = ts.tdb_jd(2458827.5)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_text(load('pck00008.tpc'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')
    assert frame.center == 301

    pt = pc.build_latlon_degrees(frame, 26.3, -46.8)
    assert pt.center == 301

    geometric = pt.at(t)
    assert geometric.center == 301

    # See the file `horizons/moon-from-moon-topos` for HORIZONS result.
    # Alas, here we achieve only about 4-5 digits of agreement.  Yes,
    # the latitude and longitude are only 4 digits each, but HORIZONS
    # considers them to be very precise values like "313.200000".  The
    # agreement here should be far greater.
    want = [1.043588965592271E-05, 3.340834944508400E-06,
            -3.848560523814720E-06]
    assert abs(geometric.position.au - want).max() < 4e-9

def test_observing_earth_from_location_on_moon():
    ts = load.timescale(builtin=True)
    t = ts.utc(2019, 12, 13)
    # t = ts.tdb_jd(2458827.5) CLEAN UP

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_text(load('pck00008.tpc'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')
    assert frame.center == 301

    pt = pc.build_latlon_degrees(frame, 26.3, -46.8)
    assert pt.center == 301

    eph = load('de421.bsp')
    astrometric = (eph['moon'] + pt).at(t).observe(eph['earth'])

    # See the file `horizons/earth-from-moon-topos` for HORIZONS result.
    print(astrometric.position.au)
    ra, dec, distance = astrometric.radec()
    ra, dec, distance = astrometric.apparent().radec()
    print(ra._degrees)
    print(dec.degrees)

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
