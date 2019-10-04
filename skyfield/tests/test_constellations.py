from skyfield.constellationlib import load_constellation_lookup

def test_constellations():
    lookup = load_constellation_lookup()

    assert lookup(24, -90) == 'Oct'
    assert lookup(0, 0) == 'Psc'
    assert lookup(4.65, 0) == 'Ori'
    assert lookup(10, 90) == 'UMi'

    assert (lookup([4.65, 10], [0, 90]) == ['Ori', 'UMi']).all()
