import numpy as np
from skyfield.api import PlanetaryConstants, T0, load
from skyfield.constants import AU_KM, AU_M
from skyfield.positionlib import ICRF

def test_frame_rotation_matrices():
    # To produce the following matrices:
    #
    # import numpy as np
    # import spiceypy as spice
    # from skyfield.constants import DAY_S
    # spice.furnsh('moon_080317.tf')
    # spice.furnsh('moon_pa_de421_1900-2050.bpc')
    # g = np.vectorize(repr)
    # print(g(spice.pxform('J2000', 'MOON_PA', (tdb - T0) * DAY_S)))
    # print(g(spice.sxform('J2000', 'MOON_PA', (tdb - T0) * DAY_S)[3:6,:3]))

    ts = load.timescale(builtin=True)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))
    frame = pc.build_frame_named('MOON_PA_DE421')
    assert frame._matrix is None  # pure segment, with no rotation applied

    # First, a moment when the angle W is nearly zero radians, so all of
    # its precision is available to the trigonometry.

    tdb = T0 - 11150

    desired_rotation = [
        [0.9994150897380264, 0.032310270603926675, 0.011203785852719871],
        [-0.034157426811763446, 0.9272642685944782, 0.37284614304233643],
        [0.0016578894811167893, -0.3730107540024127, 0.9278255379116378],
    ]
    desired_rate = [
        [-9.091699265305538e-08, 2.4682369763316187e-06,
         9.920226874449144e-07],
        [-2.660127501349402e-06, -8.620891007588608e-08,
         -2.930074158824601e-08],
        [4.245963144600506e-10, -5.067878765502927e-10,
         -2.045010122720299e-10],
    ]

    R = frame.rotation_at(ts.tdb_jd(tdb))
    assert (R == desired_rotation).all()  # Boom.

    R2, Rv = frame.rotation_and_rate_at(ts.tdb_jd(tdb))
    assert (R == R2).all()
    assert (Rv == desired_rate).all()  # Boom.

    # Second, a moment when the angle W is more than 2500 radians.

    tdb = T0
    desired_rotation = [
        [0.7840447406961362, 0.5582359944893811, 0.2713787372716964],
        [-0.6203032939745002, 0.7203957219351799, 0.31024800934393754],
        [-0.02230847532023746, -0.41158544468183367, 0.9110981032001678],
    ]
    desired_rate = [
        [-1.6512401259577911e-06, 1.9173507906460613e-06,
         8.265640603882371e-07],
        [-2.0870970217531474e-06, -1.4860137438942676e-06,
         -7.223743806558455e-07],
        [-5.817943897465853e-10, -4.4636767256698343e-10,
         -2.1589045361778893e-10],
    ]
    R = frame.rotation_at(ts.tdb_jd(tdb))
    delta = abs(R - desired_rotation)
    assert (delta < 2e-13).all()  # a few digits are lost in large W radians?

    R2, Rv = frame.rotation_and_rate_at(ts.tdb_jd(tdb))
    assert (R == R2).all()
    assert abs(Rv - desired_rate).max() < 4e-19  # About 13 digits precision.

    # Finally, a frame which is defined by a constant rotation of
    # another frame.

    tdb = T0 - 11150
    desired_rotation = [
        [0.9994268420493244, 0.03186286343877705, 0.011434392191818011],
        [-0.03382833397374688, 0.9272754001640859, 0.37284846261060917],
        [0.0012771890522186643, -0.37302156798771835, 0.927821792481783],
    ]
    desired_rate = [
        [-9.004087820606406e-08, 2.4682648578114788e-06,
         9.920321321295625e-07],
        [-2.660157295439155e-06, -8.539614878340816e-08,
         -2.8974080519384074e-08],
        [4.550211407392054e-10, -1.4469992803332584e-09,
         -5.823780954758093e-10],
    ]
    frame = pc.build_frame_named('MOON_ME_DE421')
    R = frame.rotation_at(ts.tdb_jd(tdb))
    delta = abs(R - desired_rotation)
    assert (R == desired_rotation).all()

    R2, Rv = frame.rotation_and_rate_at(ts.tdb_jd(tdb))
    assert (R == R2).all()
    assert (Rv == desired_rate).all()

def test_rotating_vector_into_frame():
    et_seconds = 259056665.1855896
    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(T0 + et_seconds / 3600. / 24.0)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    # Example from "moon_080317.tf" (the raw vector is taken from a
    # tweaked version of the script in the file that uses "J2000"):

    vector = ICRF(np.array([2.4798273371071659e+05,
                            -2.6189996683651494e+05,
                            -1.2455830876097400e+05]) / AU_KM)
    vector.t = t
    meter = 1e-3

    frame = pc.build_frame_named('MOON_PA_DE421')
    result = vector.frame_xyz(frame)
    assert max(abs(result.km - [379908.634, 33385.003, -12516.8859])) < meter

    relative_frame = pc.build_frame_named('MOON_ME_DE421')
    result = vector.frame_xyz(relative_frame)
    assert max(abs(result.km - [379892.825, 33510.118, -12661.5278])) < meter

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
    # Alas, here we achieve only about 4-5 digits of agreement.  The
    # agreement should be far greater.  Yes, the latitude and longitude
    # are only 3 digits each, but HORIZONS considers them to be very
    # precise values like "313.200000".
    want = [-1.043588965592271E-05, -3.340834944508400E-06,
            3.848560523814720E-06]
    meter = 1.0 / AU_M
    print(geometric.position.au)
    print(want)
    print(geometric.position.au - want)
    print((geometric.position.au - want) / meter)
    assert abs(geometric.position.au - want).max() < 101 * meter

    # TODO: add a "want" from the HORIZONS velocity
    #print(geometric.velocity.au_per_d)

def test_observing_earth_from_location_on_moon():
    ts = load.timescale(builtin=True)
    t = ts.utc(2019, 12, 13)

    pc = PlanetaryConstants()
    pc.read_text(load('moon_080317.tf'))
    pc.read_text(load('pck00008.tpc'))
    pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))

    frame = pc.build_frame_named('MOON_ME_DE421')
    assert frame.center == 301

    pt = pc.build_latlon_degrees(frame, 26.3, 313.2)
    assert pt.center == 301

    eph = load('de421.bsp')
    astrometric = (eph['moon'] + pt).at(t).observe(eph['earth'])

    # See the file `horizons/earth-from-moon-topos` for HORIZONS result.
    print(astrometric.position.au)
    ra, dec, distance = astrometric.radec()
    #ra, dec, distance = astrometric.apparent().radec()
    print('want: 270.2590484        -22.8079717')
    print('got:  {} {}'.format(ra._degrees, dec.degrees))
    # asdf

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
