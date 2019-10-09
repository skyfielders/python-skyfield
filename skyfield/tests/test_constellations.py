from skyfield.api import (
    Star, load, load_constellation_map, position_from_radec,
)

def test_constellations():
    ts = load.timescale(builtin=True)
    t = ts.utc(2000)
    lookup = load_constellation_map()

    assert lookup(position_from_radec(0, 0)) == 'Psc'
    assert lookup(position_from_radec(360, 90)) == 'UMi'

    # (4.65, 0) used to be in Orion, but, precession
    assert lookup(position_from_radec(4.65, 0)) == 'Eri'
    assert lookup(position_from_radec(4.75, 0.3)) == 'Ori'

    # Test vectorization.
    assert (lookup(position_from_radec([4.65, 4.75], [0, 0.3]))
            == ['Eri', 'Ori']).all()
